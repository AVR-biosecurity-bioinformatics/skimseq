#!/usr/bin/env python3

import argparse
import logging
import math
import time
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pyBigWig
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
            self.handle.write(
                f"{chrom}\t{start}\t{end}\t{self.label}\n"
            )
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

def read_bigwig_summary(path, expected_chrom_sizes=None):
    """calculate genome-wide mean coverages from BigWig header"""
    with pyBigWig.open(path) as bw:
        chrom_sizes = dict(bw.chroms())
        header = bw.header()

    if expected_chrom_sizes is not None and chrom_sizes != expected_chrom_sizes:
        raise ValueError(f"Chromosome names or lengths differ in BigWig: {path}")

    genome_size = sum(chrom_sizes.values())
    mean_coverage = float(header["sumData"]) / genome_size

    return {
        "chrom_sizes": chrom_sizes,
        "sum_data": float(header["sumData"]),
        "bases_covered": int(header.get("nBasesCovered", 0)),
        "genome_size": genome_size,
        "mean_coverage": mean_coverage,
    }

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bigwig",
        nargs="+",
        required=True,
        help="Input BigWig files.",
    )
    parser.add_argument(
        "--window-size",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--sample-batch-size",
        type=int,
        default=32,
    )
    parser.add_argument(
        "--window-batch-size",
        type=int,
        default=None,
        help=(
            "Maximum windows per FigWig call. By default, calculate this "
            "automatically from --array-memory-mib."
        ),
    )
    parser.add_argument(
        "--array-memory-mib",
        type=int,
        default=4096,
        help="Maximum size in MiB of the main dense FigWig output array.",
    )
    parser.add_argument(
        "--low",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--high",
        type=float,
        default=1.5,
    )
    parser.add_argument(
        "--prefix",
        default="cohort",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    for name in (
        "window_size",
        "threads",
        "sample_batch_size",
        "array_memory_mib",
    ):
        if getattr(args, name) <= 0:
            parser.error(
                f"--{name.replace('_', '-')} must be positive"
            )

    if (
        args.window_batch_size is not None
        and args.window_batch_size <= 0
    ):
        parser.error("--window-batch-size must be positive")

    if not (
        math.isfinite(args.low)
        and math.isfinite(args.high)
        and args.low < args.high
    ):
        parser.error(
            "Thresholds must be finite and low must be less than high"
        )

    bigwig_paths = [str(Path(path)) for path in args.bigwig]
    sample_names = [
        Path(path).name
        .removesuffix(".bigwig")
        .removesuffix(".bigWig")
        .removesuffix(".bw")
        for path in bigwig_paths
    ]

    if len(set(bigwig_paths)) != len(bigwig_paths):
        raise ValueError("Input BigWig paths must be unique")

    if len(set(sample_names)) != len(sample_names):
        raise ValueError(
            "BigWig filenames produce duplicate sample names after "
            "removing their extensions"
        )

    n_samples = len(bigwig_paths)

    summaries = [read_bigwig_summary(bigwig_paths[0])]
    chrom_sizes = summaries[0]["chrom_sizes"]

    summaries.extend(
        read_bigwig_summary(
            path,
            expected_chrom_sizes=chrom_sizes,
        )
        for path in bigwig_paths[1:]
    )

    scales = np.asarray(
        [1.0 / s["mean_coverage"] for s in summaries],
        dtype=np.float64,
    )

    for sample, path, summary, scale in zip(
        sample_names,
        bigwig_paths,
        summaries,
        scales,
    ):
        if summary["bases_covered"] != summary["genome_size"]:
            logger.warning(
                "Sample %s has %d bases represented in %s but %d bases "
                "in the reference; omitted positions are included as "
                "zero coverage when calculating the mean",
                sample,
                summary["bases_covered"],
                path,
                summary["genome_size"],
            )

    logger.info(
        "Sample %s: mean coverage %.8g; scale factor %.8g",
        sample,
        summary["mean_coverage"],
        scale,
    )

    # Retain these objects so FigWig can reuse their cached indexes.
    bigwigs = [BigWig(path) for path in bigwig_paths]

    n_windows = sum(
        (size + args.window_size - 1) // args.window_size
        for size in chrom_sizes.values()
    )

    sample_batch_size = min(
        args.sample_batch_size,
        n_samples,
    )

    # Bound the dense float32 array returned by FigWig. Leave some headroom
    # within the requested budget for temporary arrays and library overhead.
    budget_bytes = args.array_memory_mib * 1024 * 1024
    usable_budget_bytes = int(budget_bytes * 0.85)

    bytes_per_window = (
        np.dtype(np.float32).itemsize
        * sample_batch_size
        * args.window_size
    )

    if bytes_per_window > usable_budget_bytes:
        raise ValueError(
            "One window/sample batch exceeds --array-memory-mib after "
            "allowing memory headroom; reduce --sample-batch-size or "
            "increase --array-memory-mib"
        )

    memory_limited_batch_size = (
        usable_budget_bytes // bytes_per_window
    )

    # Impose a default operational cap to avoid creating extremely large
    # coordinate arrays when very small windows are requested.
    if args.window_batch_size is None:
        window_batch_size = min(
            memory_limited_batch_size,
            262_144,
        )
    else:
        window_batch_size = min(
            args.window_batch_size,
            memory_limited_batch_size,
        )

    # There is no benefit in creating batches larger than the genome.
    window_batch_size = min(
        window_batch_size,
        n_windows,
    )

    estimated_array_mib = (
        window_batch_size * bytes_per_window / 1024**2
    )
    logger.info(
        "%d samples; %d contigs; %d windows; "
        "batch limits: %d windows x %d samples; %d threads; "
        "estimated dense array: %.1f MiB",
        n_samples,
        len(chrom_sizes),
        n_windows,
        window_batch_size,
        sample_batch_size,
        args.threads,
        estimated_array_mib,
    )

    output_buffer_size = 8 * 1024 * 1024
    start_time = time.perf_counter()

    with ExitStack() as stack:
        bedgraph = stack.enter_context(
            open(
                f"{args.prefix}.normalised-coverage.bedgraph",
                "w",
                buffering=output_buffer_size,
            )
        )
        high_handle = stack.enter_context(
            open(
                f"{args.prefix}.high_depth.bed",
                "w",
                buffering=output_buffer_size,
            )
        )
        low_handle = stack.enter_context(
            open(
                f"{args.prefix}.low_depth.bed",
                "w",
                buffering=output_buffer_size,
            )
        )
        scale_handle = stack.enter_context(
            open(
                f"{args.prefix}.scale-factors.tsv",
                "w",
                buffering=output_buffer_size,
            )
        )

        scale_handle.write(
            "sample\tbigwig\tmean_coverage\tscale_factor\t"
            "sum_data\tbases_covered\tgenome_size\n"
        )

        for sample, path, summary, scale in zip(
            sample_names,
            bigwig_paths,
            summaries,
            scales,
        ):
            scale_handle.write(
                f"{sample}\t"
                f"{path}\t"
                f"{summary['mean_coverage']:.17g}\t"
                f"{scale:.17g}\t"
                f"{summary['sum_data']:.17g}\t"
                f"{summary['bases_covered']}\t"
                f"{summary['genome_size']}\n"
            )

        high_mask = MaskWriter(
            high_handle,
            "HIGH_DEPTH",
        )
        low_mask = MaskWriter(
            low_handle,
            "LOW_DEPTH",
        )

        processed = 0

        for batch_index, (chroms, starts, lengths) in enumerate(
            window_batches(
                chrom_sizes,
                args.window_size,
                window_batch_size,
            ),
            start=1,
        ):
            cohort_sum = np.zeros(
                len(starts),
                dtype=np.float64,
            )

            partial_rows = np.flatnonzero(
                lengths < args.window_size
            )

            for first in range(
                0,
                n_samples,
                sample_batch_size,
            ):
                last = min(
                    first + sample_batch_size,
                    n_samples,
                )

                values = read_windows(
                    bigwigs[first:last],
                    chroms,
                    starts,
                    width=args.window_size,
                    n_jobs=args.threads,
                    missing=0.0,
                )

                # FigWig omits the channel axis for a single input file.
                if values.ndim == 2:
                    values = values[:, np.newaxis, :]

                expected_shape = (
                    len(starts),
                    last - first,
                    args.window_size,
                )

                if values.shape != expected_shape:
                    raise ValueError(
                        f"Unexpected FigWig shape: {values.shape}; "
                        f"expected {expected_shape}"
                    )

                # Clear only padding beyond the reference end. Real non-finite
                # values within the reference remain detectable below.
                for row_index in partial_rows:
                    values[row_index, :, lengths[row_index]:] = 0.0

                window_sums = values.sum(
                    axis=2,
                    dtype=np.float64,
                )
                del values

                if not np.isfinite(window_sums).all():
                    raise ValueError(
                        f"Non-finite coverage in window batch "
                        f"{batch_index}, sample group "
                        f"{first + 1}-{last}"
                    )

                # Divide by each window's actual reference length,
                # scale each sample, then accumulate across samples.
                window_sums /= lengths[:, np.newaxis]
                cohort_sum += (
                    window_sums @ scales[first:last]
                )

            cohort_mean = cohort_sum / n_samples

            if not np.isfinite(cohort_mean).all():
                raise ValueError(
                    f"Non-finite cohort mean in window batch "
                    f"{batch_index}"
                )

            for chrom, start, length, value in zip(
                chroms,
                starts,
                lengths,
                cohort_mean,
            ):
                start = int(start)
                end = start + int(length)

                bedgraph.write(
                    f"{chrom}\t{start}\t{end}\t{value:.8g}\n"
                )
                high_mask.add(
                    chrom,
                    start,
                    end,
                    value >= args.high,
                )
                low_mask.add(
                    chrom,
                    start,
                    end,
                    value <= args.low,
                )

            processed += len(starts)

            if batch_index == 1 or batch_index % 100 == 0:
                elapsed = time.perf_counter() - start_time
                rate = processed / elapsed

                logger.info(
                    "Processed %d / %d windows at %.0f windows/s",
                    processed,
                    n_windows,
                    rate,
                )

        high_mask.flush()
        low_mask.flush()

    elapsed = time.perf_counter() - start_time

    logger.info(
        "Finished: %d windows in %.1f seconds at %.0f windows/s",
        processed,
        elapsed,
        processed / elapsed,
    )


if __name__ == "__main__":
    main()