#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
import tempfile
import numpy as np

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
        description=(
            "Call mitochondrial consensus from direct bcftools mpileup "
            "FORMAT/AD all-sites output."
        )
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

    sample_ids = [row["sample_id"] for row in rows]

    duplicates = sorted(
        sample_id
        for sample_id, count in Counter(sample_ids).items()
        if count > 1
    )

    if duplicates:
        raise ValueError(
            f"Duplicate sample IDs found in {path}: "
            f"{', '.join(duplicates)}"
        )

    return sample_ids


def split_counts(value):
    if value in {"", "."}:
        return []

    return [0 if x in {"", "."} else int(x) for x in value.split(",")]


def parse_bcftools_alleles(ref, alt):
    """
    Return alleles in the same order as bcftools FORMAT/AD:

        REF, ALT1, ALT2, ...

    The symbolic <*> allele is retained here so that allele indices remain
    aligned with the AD vector. It is removed later when parsing each sample.
    """

    ref = ref.upper()

    if alt in {"", "."}:
        return [ref]

    return [
        ref,
        *[
            allele.upper()
            for allele in alt.split(",")
            if allele
        ],
    ]


def parse_sample_field(value, bcftools_alleles):
    """
    Parse one direct bcftools mpileup FORMAT/AD sample field.

    Examples:
        460,0
        481,0,2,0
        .

    FORMAT/AD is ordered as:
        REF, ALT1, ALT2, ...

    Single-nucleotide A/C/G/T alleles contribute to SNV counts.
    Other real alleles contribute to non-SNV evidence.

    The symbolic <*> allele represents spanning-deletion evidence. It is
    excluded from the reported allele list but included in non_snv_count.

    The reported gt is a minipileup-compatible allele-index label, not a
    biological diploid genotype.
    """

    raw_counts = split_counts(value)

    if raw_counts and len(raw_counts) != len(bcftools_alleles):
        raise ValueError(
            "FORMAT/AD length does not match REF/ALT allele count: "
            f"alleles={bcftools_alleles}, AD={raw_counts}"
        )

    if not raw_counts:
        raw_counts = [0] * len(bcftools_alleles)

    counts = {base: 0 for base in BASES}
    alleles = []
    allele_counts = []
    non_snv_count = 0
    highest_supported = 0

    for allele, allele_count in zip(bcftools_alleles, raw_counts):
        allele = allele.upper()

        # <*> represents reads spanning a deletion. Include this in non-SNV
        # evidence, but exclude it from the reported allele vectors.
        if allele == "<*>":
            non_snv_count += allele_count
            continue

        output_index = len(alleles)

        alleles.append(allele)
        allele_counts.append(allele_count)

        if allele in counts:
            counts[allele] += allele_count
        elif allele not in {"", "."}:
            non_snv_count += allele_count

        if output_index > 0 and allele_count >= 1:
            highest_supported = output_index

    snv_depth = sum(counts.values())
    total_depth = snv_depth + non_snv_count

    return {
        "gt": f"0/{highest_supported}",
        "alleles": alleles,
        "allele_counts": allele_counts,
        "counts": counts,
        "snv_depth": snv_depth,
        "total_depth": total_depth,
        "non_snv_count": non_snv_count,
        "non_snv_af": (
            non_snv_count / total_depth
            if total_depth
            else 0.0
        ),
    }


def parse_counts_row(line, path, line_number, expected_cols):
    line = line.rstrip("\r\n")

    if not line:
        return None

    fields = line.split("\t")

    if len(fields) != expected_cols:
        raise ValueError(
            f"{path}:{line_number} has {len(fields)} columns, "
            f"expected {expected_cols}. Check that the sample manifest "
            f"matches the BAM order passed to bcftools mpileup."
        )

    try:
        pos = int(fields[1])
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number} has an invalid position: {fields[1]}"
        ) from exc

    return fields[0], pos, fields


def build_shifted_offset_index(
    path,
    expected_cols,
    shift_bases,
    mito_length,
    breakpoint_window,
):
    """
    Record the file offset of each shifted-pileup row needed near the
    original mitochondrial breakpoint.

    Memory complexity is O(mito_length), independent of sample count.
    """
    offsets = [None] * (mito_length + 1)

    with open(path) as handle:
        line_number = 0

        while True:
            offset = handle.tell()
            line = handle.readline()

            if not line:
                break

            line_number += 1

            parsed = parse_counts_row(
                line=line,
                path=path,
                line_number=line_number,
                expected_cols=expected_cols,
            )

            if parsed is None:
                continue

            _, shifted_pos, _ = parsed

            if not 1 <= shifted_pos <= mito_length:
                raise ValueError(
                    f"{path}:{line_number} has position {shifted_pos}, "
                    f"outside mitochondrial length {mito_length}"
                )

            original_pos = shifted_to_original_pos(
                shifted_pos=shifted_pos,
                shift_bases=shift_bases,
                mito_length=mito_length,
            )

            if is_breakpoint_pos(
                pos=original_pos,
                mito_length=mito_length,
                window=breakpoint_window,
            ):
                if offsets[original_pos] is not None:
                    raise ValueError(
                        f"{path} contains multiple records mapping to "
                        f"original position {original_pos}"
                    )

                offsets[original_pos] = offset

    return offsets


class CountsStream:
    """
    Sequential reader for a position-sorted bcftools counts file.

    Missing positions are returned as None.
    """

    def __init__(self, path, expected_cols, expected_contig=None):
        self.path = path
        self.expected_cols = expected_cols
        self.expected_contig = expected_contig
        self.handle = open(path)
        self.line_number = 0
        self.current = None
        self.previous_pos = 0
        self._advance()

    def _advance(self):
        while True:
            line = self.handle.readline()

            if not line:
                self.current = None
                return

            self.line_number += 1

            parsed = parse_counts_row(
                line=line,
                path=self.path,
                line_number=self.line_number,
                expected_cols=self.expected_cols,
            )

            if parsed is None:
                continue

            chrom, pos, fields = parsed

            if (
                self.expected_contig is not None
                and chrom != self.expected_contig
            ):
                raise ValueError(
                    f"{self.path}:{self.line_number} uses contig {chrom}, "
                    f"expected {self.expected_contig}"
                )

            if pos <= self.previous_pos:
                raise ValueError(
                    f"{self.path}:{self.line_number} is not strictly "
                    f"position-sorted: {pos} follows {self.previous_pos}"
                )

            self.previous_pos = pos
            self.current = chrom, pos, fields
            return

    def get(self, pos):
        while self.current is not None and self.current[1] < pos:
            self._advance()

        if self.current is not None and self.current[1] == pos:
            result = self.current
            self._advance()
            return result

        return None

    def close(self):
        self.handle.close()


def read_shifted_row(handle, offset, path, expected_cols):
    if offset is None:
        return None

    handle.seek(offset)
    line = handle.readline()

    return parse_counts_row(
        line=line,
        path=path,
        line_number=0,
        expected_cols=expected_cols,
    )


def shifted_to_original_pos(shifted_pos, shift_bases, mito_length):
    return ((shifted_pos + shift_bases - 1) % mito_length) + 1



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
    n_samples = len(samples)
    expected_cols = 4 + n_samples

    shifted_offsets = build_shifted_offset_index(
        path=args.shifted_counts,
        expected_cols=expected_cols,
        shift_bases=args.shift_bases,
        mito_length=mito_length,
        breakpoint_window=args.breakpoint_window,
    )

    call_fields = [
        "sample_id",
        "contig",
        "pos",
        "ref",
        "alleles",
        "gt",
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

    # Compact per-sample state retained in memory.
    consensus = [bytearray(b"N" * mito_length) for _ in samples]
    total_depth_sum = [0] * n_samples
    snv_depth_sum = [0] * n_samples
    covered_bases = [0] * n_samples
    any_non_snv_sites = [0] * n_samples
    shifted_source_bases = [0] * n_samples
    filters = [Counter() for _ in samples]

    original_stream = CountsStream(
        path=args.original_counts,
        expected_cols=expected_cols,
        expected_contig=contig,
    )

    try:
        with tempfile.TemporaryDirectory(
            prefix="mito_consensus_"
        ) as temp_dir, open(
            args.shifted_counts
        ) as shifted_handle, open(
            args.out_calls,
            "w",
            newline="",
        ) as calls_out:

            # Disk-backed depth matrices for exact median calculation.
            # Rows are positions and columns are samples.
            total_depths = np.memmap(
                f"{temp_dir}/total_depths.uint32",
                dtype=np.uint32,
                mode="w+",
                shape=(mito_length, n_samples),
            )

            snv_depths = np.memmap(
                f"{temp_dir}/snv_depths.uint32",
                dtype=np.uint32,
                mode="w+",
                shape=(mito_length, n_samples),
            )

            calls_writer = csv.DictWriter(
                calls_out,
                delimiter="\t",
                fieldnames=call_fields,
            )
            calls_writer.writeheader()

            for pos in range(1, mito_length + 1):
                original_row = original_stream.get(pos)

                row = original_row
                source = "original"
                shifted_pos = "."

                # Use the shifted pileup around the circular breakpoint when
                # a corresponding shifted observation is available.
                if is_breakpoint_pos(
                    pos=pos,
                    mito_length=mito_length,
                    window=args.breakpoint_window,
                ):
                    shifted_row = read_shifted_row(
                        handle=shifted_handle,
                        offset=shifted_offsets[pos],
                        path=args.shifted_counts,
                        expected_cols=expected_cols,
                    )

                    if shifted_row is not None:
                        row = shifted_row
                        source = "shifted"
                        shifted_pos = shifted_row[1]

                if row is None:
                    ref = refseq[pos - 1]
                    bcftools_alleles = None
                    sample_values = None
                    source = "missing"
                    output_ref = ref if ref in BASES else "N"

                else:
                    _, row_pos, fields = row

                    ref = fields[2].upper()
                    alt = fields[3]

                    bcftools_alleles = parse_bcftools_alleles(
                        ref=ref,
                        alt=alt,
                    )

                    sample_values = fields[4:]
                    output_ref = ref if ref in BASES else "N"

                    if source == "original" and row_pos != pos:
                        raise ValueError(
                            f"Original pileup position mismatch: "
                            f"requested {pos}, found {row_pos}"
                        )

                for sample_index, sample_id in enumerate(samples):
                    if sample_values is None:
                        obs = empty_observation(
                            contig=contig,
                            pos=pos,
                            ref_base=refseq[pos - 1],
                        )

                    else:
                        parsed = parse_sample_field(
                            value=sample_values[sample_index],
                            bcftools_alleles=bcftools_alleles,
                        )

                        obs = {
                            "chrom": contig,
                            "pos": pos,
                            "ref": output_ref,
                            **parsed,
                        }

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

                    total_depth = called["total_depth"]
                    snv_depth = called["snv_depth"]

                    if total_depth > np.iinfo(np.uint32).max:
                        raise OverflowError(
                            f"Total depth exceeds uint32 at "
                            f"{sample_id}:{pos}: {total_depth}"
                        )

                    if snv_depth > np.iinfo(np.uint32).max:
                        raise OverflowError(
                            f"SNV depth exceeds uint32 at "
                            f"{sample_id}:{pos}: {snv_depth}"
                        )

                    consensus[sample_index][pos - 1] = ord(called["call"])

                    total_depths[pos - 1, sample_index] = total_depth
                    snv_depths[pos - 1, sample_index] = snv_depth

                    total_depth_sum[sample_index] += total_depth
                    snv_depth_sum[sample_index] += snv_depth

                    if total_depth >= args.min_depth:
                        covered_bases[sample_index] += 1

                    filters[sample_index][called["filter"]] += 1

                    if obs["non_snv_count"] > 0:
                        any_non_snv_sites[sample_index] += 1

                    if source == "shifted":
                        shifted_source_bases[sample_index] += 1

                    calls_writer.writerow({
                        "sample_id": sample_id,
                        "contig": contig,
                        "pos": pos,
                        "ref": obs["ref"],
                        "alleles": ",".join(obs["alleles"]),
                        "gt": obs["gt"],
                        "allele_counts": ",".join(
                            map(str, obs["allele_counts"])
                        ),
                        "a_count": called["a_count"],
                        "c_count": called["c_count"],
                        "g_count": called["g_count"],
                        "t_count": called["t_count"],
                        "snv_depth": snv_depth,
                        "total_depth": total_depth,
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
                        "shifted_pos": (
                            shifted_pos if source == "shifted" else "."
                        ),
                    })

            total_depths.flush()
            snv_depths.flush()

            with open(
                args.out_fasta,
                "w",
            ) as fasta_out, open(
                args.out_qc,
                "w",
                newline="",
            ) as qc_out:

                qc_writer = csv.DictWriter(
                    qc_out,
                    delimiter="\t",
                    fieldnames=qc_fields,
                )
                qc_writer.writeheader()

                for sample_index, sample_id in enumerate(samples):
                    seq = consensus[sample_index].decode("ascii")
                    sample_filters = filters[sample_index]

                    mean_total_depth = (
                        total_depth_sum[sample_index] / mito_length
                    )
                    mean_snv_depth = (
                        snv_depth_sum[sample_index] / mito_length
                    )

                    median_total_depth = float(
                        np.median(total_depths[:, sample_index])
                    )
                    median_snv_depth = float(
                        np.median(snv_depths[:, sample_index])
                    )

                    n_bases = seq.count("N")

                    fasta_out.write(
                        f">{sample_id} {contig}:1-{mito_length}\n"
                    )
                    fasta_out.write(wrap_fasta(seq) + "\n")

                    qc_writer.writerow({
                        "sample_id": sample_id,
                        "mito_length": mito_length,
                        "mean_total_depth": f"{mean_total_depth:.3f}",
                        "median_total_depth": (
                            f"{median_total_depth:.3f}"
                        ),
                        "mean_snv_depth": f"{mean_snv_depth:.3f}",
                        "median_snv_depth": (
                            f"{median_snv_depth:.3f}"
                        ),
                        "covered_bases": covered_bases[sample_index],
                        "covered_fraction": (
                            f"{covered_bases[sample_index] / mito_length:.6f}"
                        ),
                        "pass_bases": sample_filters["PASS"],
                        "n_bases": n_bases,
                        "n_fraction": f"{n_bases / mito_length:.6f}",
                        "low_depth_sites": (
                            sample_filters["LOW_DEPTH"]
                        ),
                        "low_major_af_sites": (
                            sample_filters["LOW_MAJOR_AF"]
                        ),
                        "mixed_iupac_sites": (
                            sample_filters["MIXED_IUPAC"]
                        ),
                        "no_base_support_sites": (
                            sample_filters["NO_BASE_SUPPORT"]
                        ),
                        "non_snv_evidence_sites": (
                            sample_filters["NON_SNV_EVIDENCE"]
                        ),
                        "any_non_snv_evidence_sites": (
                            any_non_snv_sites[sample_index]
                        ),
                        "shifted_source_bases": (
                            shifted_source_bases[sample_index]
                        ),
                    })

            # Explicitly release memory-map handles before TemporaryDirectory
            # attempts to remove the backing files.
            del total_depths
            del snv_depths

    finally:
        original_stream.close()


if __name__ == "__main__":
    main()