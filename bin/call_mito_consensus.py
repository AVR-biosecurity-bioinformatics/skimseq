#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict

from pyfaidx import Fasta

BASES = ("A", "C", "G", "T")

IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
    frozenset(("A", "C", "G")): "V",
    frozenset(("A", "C", "T")): "H",
    frozenset(("A", "G", "T")): "D",
    frozenset(("C", "G", "T")): "B",
    frozenset(("A", "C", "G", "T")): "N",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Call mitochondrial consensus from minipileup -C all-sites output."
    )

    p.add_argument("--samples", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--original-counts", required=True)
    p.add_argument("--shifted-counts", required=True)

    p.add_argument("--shift-bases", type=int, required=True)
    p.add_argument("--breakpoint-window", type=int, default=500)

    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--major-af", type=float, default=0.80)
    p.add_argument("--mixed-min-af", type=float, default=0.20)
    p.add_argument("--min-minor-depth", type=int, default=3)
    p.add_argument(
        "--max-non-snv-af",
        type=float,
        default=0.20,
        help=(
            "Maximum allowed fraction of non-SNV evidence. Sites with "
            "non_snv_count / total_depth above this threshold are masked as N."
        ),
    )

    p.add_argument(
        "--het-mode",
        choices=("N", "iupac"),
        default="N",
        help="How to handle mixed SNV sites below --major-af.",
    )

    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-calls", required=True)
    p.add_argument("--out-qc", required=True)

    return p.parse_args()


def read_reference(path):
    fa = Fasta(path, as_raw=True, sequence_always_upper=True)
    contigs = list(fa.keys())

    if len(contigs) != 1:
        raise ValueError(
            f"Expected one mitochondrial FASTA record, found {len(contigs)}: {contigs}"
        )

    contig = contigs[0]
    seq = str(fa[contig]).upper()

    return contig, seq


def read_samples(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"input_index", "sample_id"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )

        rows = sorted(reader, key=lambda r: int(r["input_index"]))

    if not rows:
        raise ValueError(f"No samples found in {path}")

    return [row["sample_id"] for row in rows]


def split_counts(value):
    if value in {"", "."}:
        return []

    return [0 if x in {"", "."} else int(x) for x in value.split(",")]


def parse_alleles(value):
    alleles = [x.upper() for x in value.split(",") if x]

    if not alleles:
        raise ValueError("Empty allele column")

    return alleles


def parse_sample_field(value, alleles):
    """
    Parse minipileup -C sample field.

    Expected format:
        GT:FWD_COUNTS:REV_COUNTS

    Examples:
        0/0:12:7
        0/1:12,3:7,2
        ./.:0:0

    Alleles that are exactly A/C/G/T contribute to SNV/base counts.
    All other allele strings contribute to non_snv_count.
    """

    parts = value.split(":")

    gt = parts[0] if len(parts) > 0 else "./."
    fwd = split_counts(parts[1]) if len(parts) > 1 else []
    rev = split_counts(parts[2]) if len(parts) > 2 else []

    counts = {base: 0 for base in BASES}
    allele_counts = []
    non_snv_count = 0

    n = max(len(alleles), len(fwd), len(rev))

    for i in range(n):
        allele = alleles[i].upper() if i < len(alleles) else "."
        fwd_count = fwd[i] if i < len(fwd) else 0
        rev_count = rev[i] if i < len(rev) else 0
        total = fwd_count + rev_count

        allele_counts.append(total)

        if allele in counts:
            counts[allele] += total
        else:
            non_snv_count += total

    snv_depth = sum(counts.values())
    total_depth = snv_depth + non_snv_count

    return {
        "gt": gt,
        "fwd_raw": ",".join(map(str, fwd)) if fwd else ".",
        "rev_raw": ",".join(map(str, rev)) if rev else ".",
        "allele_counts": allele_counts,
        "counts": counts,
        "snv_depth": snv_depth,
        "total_depth": total_depth,
        "non_snv_count": non_snv_count,
        "non_snv_af": non_snv_count / total_depth if total_depth else 0.0,
    }


def read_minipileup(path, samples):
    """
    Read minipileup all-sites table:

        chrom pos ref alleles sample1 sample2 ...

    No header is expected.
    """

    expected_cols = 4 + len(samples)
    data = defaultdict(lambda: defaultdict(dict))

    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != expected_cols:
                raise ValueError(
                    f"{path}:{line_number} has {len(fields)} columns, "
                    f"expected {expected_cols}. "
                    f"Check that all.samples.tsv matches the BAM order passed to minipileup."
                )

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2].upper()
            alleles = parse_alleles(fields[3])

            if ref not in BASES:
                ref = "N"

            for sample_id, sample_value in zip(samples, fields[4:]):
                parsed = parse_sample_field(sample_value, alleles)

                data[sample_id][chrom][pos] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alleles": alleles,
                    **parsed,
                }

    return data


def shifted_to_original_pos(shifted_pos, shift_bases, mito_length):
    return ((shifted_pos + shift_bases - 1) % mito_length) + 1


def remap_shifted_data(shifted_data, original_contig, mito_length, shift_bases):
    remapped = defaultdict(lambda: defaultdict(dict))

    for sample_id, by_contig in shifted_data.items():
        if len(by_contig) != 1:
            raise ValueError(
                f"Expected one shifted mitochondrial contig for {sample_id}, "
                f"found {len(by_contig)}"
            )

        shifted_contig = next(iter(by_contig))

        for shifted_pos, obs in by_contig[shifted_contig].items():
            original_pos = shifted_to_original_pos(
                shifted_pos=shifted_pos,
                shift_bases=shift_bases,
                mito_length=mito_length,
            )

            new_obs = dict(obs)
            new_obs["chrom"] = original_contig
            new_obs["pos"] = original_pos
            new_obs["shifted_chrom"] = shifted_contig
            new_obs["shifted_pos"] = shifted_pos

            remapped[sample_id][original_contig][original_pos] = new_obs

    return remapped


def is_breakpoint_pos(pos, mito_length, window):
    return window > 0 and (pos <= window or pos > mito_length - window)


def call_major_allele(
    counts,
    non_snv_count,
    min_depth,
    major_af,
    mixed_min_af,
    min_minor_depth,
    max_non_snv_af,
    het_mode,
):
    snv_depth = sum(counts.values())
    total_depth = snv_depth + non_snv_count
    non_snv_af = non_snv_count / total_depth if total_depth else 0.0

    ranked = sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0]),
    )

    major_base, major_count = ranked[0]
    second_base, second_count = ranked[1]

    major_fraction = major_count / snv_depth if snv_depth else 0.0
    second_fraction = second_count / snv_depth if snv_depth else 0.0

    if total_depth < min_depth:
        call = "N"
        filt = "LOW_DEPTH"

    elif non_snv_af > max_non_snv_af:
        call = "N"
        filt = "NON_SNV_EVIDENCE"

    elif snv_depth == 0 or major_count == 0:
        call = "N"
        filt = "NO_BASE_SUPPORT"

    elif major_fraction >= major_af:
        call = major_base
        filt = "PASS"

    elif (
        het_mode == "iupac"
        and second_count >= min_minor_depth
        and second_fraction >= mixed_min_af
    ):
        call = IUPAC.get(frozenset((major_base, second_base)), "N")
        filt = "MIXED_IUPAC"

    else:
        call = "N"
        filt = "LOW_MAJOR_AF"

    return {
        "call": call,
        "filter": filt,
        "total_depth": total_depth,
        "snv_depth": snv_depth,
        "non_snv_count": non_snv_count,
        "non_snv_af": non_snv_af,
        "a_count": counts["A"],
        "c_count": counts["C"],
        "g_count": counts["G"],
        "t_count": counts["T"],
        "major_base": major_base,
        "major_count": major_count,
        "major_af": major_fraction,
        "second_base": second_base,
        "second_count": second_count,
        "second_af": second_fraction,
    }


def wrap_fasta(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def empty_observation(contig, pos, ref_base):
    if ref_base not in BASES:
        ref_base = "N"

    return {
        "chrom": contig,
        "pos": pos,
        "ref": ref_base,
        "alleles": [ref_base],
        "gt": "./.",
        "fwd_raw": ".",
        "rev_raw": ".",
        "allele_counts": [0],
        "counts": {base: 0 for base in BASES},
        "snv_depth": 0,
        "total_depth": 0,
        "non_snv_count": 0,
        "non_snv_af": 0.0,
    }


def main():
    args = parse_args()

    contig, refseq = read_reference(args.reference)
    mito_length = len(refseq)

    if not 0 <= args.shift_bases < mito_length:
        raise ValueError(
            f"--shift-bases must be >= 0 and < mitochondrial length "
            f"({mito_length}); got {args.shift_bases}"
        )

    if args.breakpoint_window < 0:
        raise ValueError("--breakpoint-window must be >= 0")

    if not 0 <= args.major_af <= 1:
        raise ValueError("--major-af must be between 0 and 1")

    if not 0 <= args.mixed_min_af <= 1:
        raise ValueError("--mixed-min-af must be between 0 and 1")

    if not 0 <= args.max_non_snv_af <= 1:
        raise ValueError("--max-non-snv-af must be between 0 and 1")

    samples = read_samples(args.samples)

    original = read_minipileup(args.original_counts, samples)
    shifted_raw = read_minipileup(args.shifted_counts, samples)
    shifted = remap_shifted_data(
        shifted_data=shifted_raw,
        original_contig=contig,
        mito_length=mito_length,
        shift_bases=args.shift_bases,
    )

    call_fields = [
        "sample_id",
        "contig",
        "pos",
        "ref",
        "alleles",
        "gt",
        "fwd_counts",
        "rev_counts",
        "allele_counts",
        "a_count",
        "c_count",
        "g_count",
        "t_count",
        "snv_depth",
        "total_depth",
        "non_snv_count",
        "non_snv_af",
        "major_base",
        "major_count",
        "major_af",
        "second_base",
        "second_count",
        "second_af",
        "call",
        "filter",
        "source_pileup",
        "shifted_pos",
    ]

    qc_fields = [
        "sample_id",
        "mito_length",
        "mean_total_depth",
        "median_total_depth",
        "mean_snv_depth",
        "median_snv_depth",
        "covered_bases",
        "covered_fraction",
        "pass_bases",
        "n_bases",
        "n_fraction",
        "low_depth_sites",
        "low_major_af_sites",
        "mixed_iupac_sites",
        "no_base_support_sites",
        "non_snv_evidence_sites",
        "any_non_snv_evidence_sites",
        "shifted_source_bases",
    ]

    with open(args.out_fasta, "w") as fasta_out, \
         open(args.out_calls, "w", newline="") as calls_out, \
         open(args.out_qc, "w", newline="") as qc_out:

        calls_writer = csv.DictWriter(calls_out, delimiter="\t", fieldnames=call_fields)
        qc_writer = csv.DictWriter(qc_out, delimiter="\t", fieldnames=qc_fields)

        calls_writer.writeheader()
        qc_writer.writeheader()

        for sample_id in samples:
            consensus = []
            total_depths = []
            snv_depths = []
            filters = Counter()
            shifted_source_bases = 0
            any_non_snv_evidence_sites = 0

            for pos in range(1, mito_length + 1):
                source = "original"
                shifted_pos = "."

                obs = original.get(sample_id, {}).get(contig, {}).get(pos)

                if is_breakpoint_pos(pos, mito_length, args.breakpoint_window):
                    shifted_obs = shifted.get(sample_id, {}).get(contig, {}).get(pos)
                    if shifted_obs is not None:
                        obs = shifted_obs
                        source = "shifted"
                        shifted_pos = shifted_obs.get("shifted_pos", ".")
                        shifted_source_bases += 1

                if obs is None:
                    obs = empty_observation(contig, pos, refseq[pos - 1])
                    source = "missing"

                called = call_major_allele(
                    counts=obs["counts"],
                    non_snv_count=obs["non_snv_count"],
                    min_depth=args.min_depth,
                    major_af=args.major_af,
                    mixed_min_af=args.mixed_min_af,
                    min_minor_depth=args.min_minor_depth,
                    max_non_snv_af=args.max_non_snv_af,
                    het_mode=args.het_mode,
                )

                consensus.append(called["call"])
                total_depths.append(called["total_depth"])
                snv_depths.append(called["snv_depth"])
                filters[called["filter"]] += 1

                if obs["non_snv_count"] > 0:
                    any_non_snv_evidence_sites += 1

                calls_writer.writerow({
                    "sample_id": sample_id,
                    "contig": contig,
                    "pos": pos,
                    "ref": obs["ref"],
                    "alleles": ",".join(obs["alleles"]),
                    "gt": obs["gt"],
                    "fwd_counts": obs["fwd_raw"],
                    "rev_counts": obs["rev_raw"],
                    "allele_counts": ",".join(map(str, obs["allele_counts"])),
                    "a_count": called["a_count"],
                    "c_count": called["c_count"],
                    "g_count": called["g_count"],
                    "t_count": called["t_count"],
                    "snv_depth": called["snv_depth"],
                    "total_depth": called["total_depth"],
                    "non_snv_count": called["non_snv_count"],
                    "non_snv_af": f"{called['non_snv_af']:.6f}",
                    "major_base": called["major_base"],
                    "major_count": called["major_count"],
                    "major_af": f"{called['major_af']:.6f}",
                    "second_base": called["second_base"],
                    "second_count": called["second_count"],
                    "second_af": f"{called['second_af']:.6f}",
                    "call": called["call"],
                    "filter": called["filter"],
                    "source_pileup": source,
                    "shifted_pos": shifted_pos,
                })

            seq = "".join(consensus)

            sorted_total_depths = sorted(total_depths)
            sorted_snv_depths = sorted(snv_depths)

            mean_total_depth = (
                sum(total_depths) / len(total_depths)
                if total_depths else 0.0
            )
            mean_snv_depth = (
                sum(snv_depths) / len(snv_depths)
                if snv_depths else 0.0
            )

            if not sorted_total_depths:
                median_total_depth = 0.0
            elif len(sorted_total_depths) % 2:
                median_total_depth = sorted_total_depths[len(sorted_total_depths) // 2]
            else:
                i = len(sorted_total_depths) // 2
                median_total_depth = (
                    sorted_total_depths[i - 1] + sorted_total_depths[i]
                ) / 2

            if not sorted_snv_depths:
                median_snv_depth = 0.0
            elif len(sorted_snv_depths) % 2:
                median_snv_depth = sorted_snv_depths[len(sorted_snv_depths) // 2]
            else:
                i = len(sorted_snv_depths) // 2
                median_snv_depth = (
                    sorted_snv_depths[i - 1] + sorted_snv_depths[i]
                ) / 2

            covered_bases = sum(d >= args.min_depth for d in total_depths)
            n_bases = seq.count("N")

            fasta_out.write(f">{sample_id} {contig}:1-{mito_length}\n")
            fasta_out.write(wrap_fasta(seq) + "\n")

            qc_writer.writerow({
                "sample_id": sample_id,
                "mito_length": mito_length,
                "mean_total_depth": f"{mean_total_depth:.3f}",
                "median_total_depth": f"{median_total_depth:.3f}",
                "mean_snv_depth": f"{mean_snv_depth:.3f}",
                "median_snv_depth": f"{median_snv_depth:.3f}",
                "covered_bases": covered_bases,
                "covered_fraction": f"{covered_bases / mito_length:.6f}",
                "pass_bases": filters["PASS"],
                "n_bases": n_bases,
                "n_fraction": f"{n_bases / mito_length:.6f}",
                "low_depth_sites": filters["LOW_DEPTH"],
                "low_major_af_sites": filters["LOW_MAJOR_AF"],
                "mixed_iupac_sites": filters["MIXED_IUPAC"],
                "no_base_support_sites": filters["NO_BASE_SUPPORT"],
                "non_snv_evidence_sites": filters["NON_SNV_EVIDENCE"],
                "any_non_snv_evidence_sites": any_non_snv_evidence_sites,
                "shifted_source_bases": shifted_source_bases,
            })


if __name__ == "__main__":
    main()