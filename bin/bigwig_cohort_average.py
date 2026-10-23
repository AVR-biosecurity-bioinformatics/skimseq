#!/usr/bin/env python3

import argparse
import csv
import logging
import math
from contextlib import ExitStack

import numpy as np
from figwig import BigWigReader as BigWig
from figwig import read_bigwig as read_windows

logger = logging.getLogger(__name__)
class MaskWriter:
    """Merge adjacent qualifying windows while writing a labelled BED."""

    def __init__(self, handle, label):
        self.handle = handle
        self.label = label
        self.pending = None

    def flush(self):
        if self.pending is not None:
            chrom, start, end = self.pending
            self.handle.write(f"{chrom}\t{start}\t{end}\t{self.label}\n")
            self.pending = None

    def add(self, chrom, start, end, keep):
        if not keep:
            self.flush()
        elif (
            self.pending is not None
            and self.pending[0] == chrom
            and self.pending[2] == start
        ):
            self.pending = (chrom, self.pending[1], end)
        else:
            self.flush()
            self.pending = (chrom, start, end)


def window_batches(chrom_sizes, width, batch_size):
    """Yield bounded batches without constructing a genome-wide window table."""

    chroms, starts, lengths = [], [], []

    for chrom, size in chrom_sizes.items():
        for start in range(0, size, width):
            chroms.append(chrom)
            starts.append(start)
            lengths.append(min(width, size - start))

            if len(starts) == batch_size:
                yield (
                    np.asarray(chroms),
                    np.asarray(starts, dtype=np.int64),
                    np.asarray(lengths, dtype=np.int64),
                )
                chroms, starts, lengths = [], [], []

    if starts:
        yield (
            np.asarray(chroms),
            np.asarray(starts, dtype=np.int64),
            np.asarray(lengths, dtype=np.int64),
        )


def read_manifest(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "bigwig", "mean_coverage", "scale_factor"}

        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Manifest must contain columns: {sorted(required)}")

        rows = list(reader)

    if not rows:
        raise ValueError("Manifest contains no samples")

    samples = [row["sample"] for row in rows]
    if len(set(samples)) != len(samples):
        raise ValueError("Sample names must be unique")

    for row in rows:
        mean = float(row["mean_coverage"])
        scale = float(row["scale_factor"])

        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(f"Invalid mean coverage for {row['sample']}: {mean}")

        if not math.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid scale factor for {row['sample']}: {scale}")

        if not math.isclose(scale, 1.0 / mean, rel_tol=1e-12):
            raise ValueError(f"Inconsistent scale factor for {row['sample']}")

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--window-size", required=True, type=int)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sample-batch-size", type=int, default=8)
    parser.add_argument("--window-batch-size", type=int, default=4096)
    parser.add_argument("--array-memory-mib", type=int, default=256)
    parser.add_argument("--low", type=float, default=0.5)
    parser.add_argument("--high", type=float, default=1.5)
    parser.add_argument("--prefix", default="cohort")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    for name in (
        "window_size",
        "threads",
        "sample_batch_size",
        "window_batch_size",
        "array_memory_mib",
    ):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")

    if not (
        math.isfinite(args.low)
        and math.isfinite(args.high)
        and args.low < args.high
    ):
        parser.error("Thresholds must be finite and low must be less than high")

    rows = read_manifest(args.manifest)
    n_samples = len(rows)
    scales = np.asarray(
        [float(row["scale_factor"]) for row in rows],
        dtype=np.float64,
    )

    # Retain these objects so figwig can reuse their cached indexes.
    bigwigs = [BigWig(row["bigwig"]) for row in rows]
    chrom_sizes = dict(bigwigs[0].chrom_sizes)

    if not chrom_sizes or any(int(size) <= 0 for size in chrom_sizes.values()):
        raise ValueError("The first BigWig has an empty or invalid chromosome header")

    # Be strict rather than silently treating missing chromosomes as zero.
    for row, bw in zip(rows[1:], bigwigs[1:]):
        if dict(bw.chrom_sizes) != chrom_sizes:
            raise ValueError(
                f"Chromosome names/lengths differ for sample {row['sample']}"
            )

    sample_batch_size = min(args.sample_batch_size, n_samples)

    # Bound the main dense float32 read array, not total process memory.
    budget_bytes = args.array_memory_mib * 1024 * 1024
    bytes_per_window = (
        np.dtype(np.float32).itemsize
        * sample_batch_size
        * args.window_size
    )

    if bytes_per_window > budget_bytes:
        raise ValueError(
            "One window/sample batch exceeds --array-memory-mib; "
            "reduce --sample-batch-size or increase the memory budget"
        )

    window_batch_size = min(
        args.window_batch_size,
        budget_bytes // bytes_per_window,
    )

    n_windows = sum(
        (size + args.window_size - 1) // args.window_size
        for size in chrom_sizes.values()
    )

    logging.info(
        "%d samples; %d contigs; %d windows; "
        "batch limits: %d windows x %d samples; %d threads",
        n_samples,
        len(chrom_sizes),
        n_windows,
        window_batch_size,
        sample_batch_size,
        args.threads,
    )

    with ExitStack() as stack:
        bedgraph = stack.enter_context(
            open(f"{args.prefix}.normalised-coverage.bedgraph", "w")
        )
        high_handle = stack.enter_context(
            open(f"{args.prefix}.high_depth.bed", "w")
        )
        low_handle = stack.enter_context(
            open(f"{args.prefix}.low_depth.bed", "w")
        )

        high_mask = MaskWriter(high_handle, "HIGH_DEPTH")
        low_mask = MaskWriter(low_handle, "LOW_DEPTH")
        processed = 0

        for batch_index, (chroms, starts, lengths) in enumerate(
            window_batches(
                chrom_sizes,
                args.window_size,
                window_batch_size,
            ),
            start=1,
        ):
            cohort_sum = np.zeros(len(starts), dtype=np.float64)
            partial_rows = np.flatnonzero(lengths < args.window_size)

            for first in range(0, n_samples, sample_batch_size):
                last = min(first + sample_batch_size, n_samples)

                values = read_windows(
                    bigwigs[first:last],
                    chroms,
                    starts,
                    width=args.window_size,
                    n_jobs=args.threads,
                    missing=0.0,
                )

                # figwig omits the channel axis for a single input file.
                if values.ndim == 2:
                    values = values[:, np.newaxis, :]

                expected_shape = (
                    len(starts),
                    last - first,
                    args.window_size,
                )
                if values.shape != expected_shape:
                    raise ValueError(
                        f"Unexpected figwig shape: {values.shape}; "
                        f"expected {expected_shape}"
                    )

                # Clear only out-of-reference padding. Real non-finite
                # values within the reference remain detectable below.
                for row_index in partial_rows:
                    values[row_index, :, int(lengths[row_index]):] = 0.0

                window_sums = values.sum(axis=2, dtype=np.float64)
                del values

                if not np.isfinite(window_sums).all():
                    raise ValueError(
                        f"Non-finite coverage in window batch {batch_index}, "
                        f"sample group {first + 1}-{last}"
                    )

                # Divide by each window's actual reference length,
                # scale each sample, then accumulate across samples.
                window_sums /= lengths[:, np.newaxis]
                cohort_sum += window_sums @ scales[first:last]

            cohort_mean = cohort_sum / n_samples

            if not np.isfinite(cohort_mean).all():
                raise ValueError(
                    f"Non-finite cohort mean in window batch {batch_index}"
                )

            for chrom, start, length, value in zip(
                chroms, starts, lengths, cohort_mean
            ):
                start = int(start)
                end = start + int(length)

                bedgraph.write(f"{chrom}\t{start}\t{end}\t{value:.17g}\n")
                high_mask.add(chrom, start, end, value >= args.high)
                low_mask.add(chrom, start, end, value <= args.low)

            processed += len(starts)
            if batch_index == 1 or batch_index % 100 == 0:
                logging.info("Processed %d / %d windows", processed, n_windows)

        high_mask.flush()
        low_mask.flush()

    logging.info("Finished: %d windows", processed)


if __name__ == "__main__":
    main()