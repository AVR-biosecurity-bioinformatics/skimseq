#!/usr/bin/env python3

import argparse
import logging
import math
import time
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
    """Calculate genome-wide mean coverage from a BigWig header."""

    with pyBigWig.open(path) as bw:
        chrom_sizes = dict(bw.chroms())
        header = bw.header()

    if not chrom_sizes or any(size <= 0 for size in chrom_sizes.values()):
        raise ValueError(f"BigWig has invalid chromosome sizes: {path}")

    if expected_chrom_sizes is not None and chrom_sizes != expected_chrom_sizes:
        raise ValueError(
            f"Chromosome names or lengths differ in BigWig: {path}"
        )

    if not header or "sumData" not in header:
        raise ValueError(f"BigWig has no stored sumData summary: {path}")

    genome_size = sum(chrom_sizes.values())
    sum_data = float(header["sumData"])
    mean_coverage = sum_data / genome_size

    return {
        "chrom_sizes": chrom_sizes,
        "mean_coverage": mean_coverage,
    }

def write_depth_masks(
        chrom_sizes,
        width,
        values,
        low_path,
        high_path,
        low,
        high,
    ):
    """Create merged low- and high-depth masks from cohort window means."""

    with (
        open(low_path, "w") as low_handle,
        open(high_path, "w") as high_handle,
    ):
        low_mask = MaskWriter(low_handle, "LOW_DEPTH")
        high_mask = MaskWriter(high_handle, "HIGH_DEPTH")

        value_index = 0

        for chrom, size in chrom_sizes.items():
            for start in range(0, size, width):
                end = min(start + width, size)
                value = values[value_index]
                value_index += 1

                low_mask.add(chrom, start, end, value <= low)
                high_mask.add(chrom, start, end, value >= high)

        low_mask.flush()
        high_mask.flush()

    if value_index != len(values):
        raise ValueError(
            f"Generated {value_index} windows but received "
            f"{len(values)} cohort means"
        )

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
        "--low-percentile",
        type=float,
        default=1.0,
        help="Mask windows at or below this percentile.",
    )

    parser.add_argument(
        "--high-percentile",
        type=float,
        default=99.0,
        help="Mask windows at or above this percentile.",
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
        math.isfinite(args.low_percentile)
        and 0.0 <= args.low_percentile < 100.0
    ):
        parser.error("--low-percentile must be between 0 and 100")

    if not (
        math.isfinite(args.high_percentile)
        and 0.0 < args.high_percentile <= 100.0
    ):
        parser.error("--high-percentile must be between 0 and 100")

    if args.low_percentile >= args.high_percentile:
        parser.error(
            "--low-percentile must be less than --high-percentile"
        )

    bigwig_paths = args.bigwig

    if len(set(bigwig_paths)) != len(bigwig_paths):
        raise ValueError("Input BigWig paths must be unique")

    # Read BigWig summaries and require a consistent reference.
    summaries = []
    chrom_sizes = None

    for path in bigwig_paths:
        summary = read_bigwig_summary(
            path,
            expected_chrom_sizes=chrom_sizes,
        )

        if chrom_sizes is None:
            chrom_sizes = summary["chrom_sizes"]

        summaries.append(summary)

    # BigWigs with non-positive or non-finite mean coverage cannot be normalised.
    valid = []

    for path, summary in zip(bigwig_paths, summaries):
        mean_coverage = summary["mean_coverage"]

        if not math.isfinite(mean_coverage) or mean_coverage <= 0.0:
            logger.warning(
                "Excluding %s because it has invalid mean coverage: %.8g",
                Path(path).name,
                mean_coverage,
            )
            continue

        valid.append((path, summary))

    if not valid:
        raise ValueError(
            "No input BigWigs have finite, positive mean coverage"
        )

    bigwig_paths = [path for path, _ in valid]
    summaries = [summary for _, summary in valid]

    n_samples = len(bigwig_paths)

    mean_coverages = np.asarray(
        [summary["mean_coverage"] for summary in summaries],
        dtype=np.float64,
    )
    scales = 1.0 / mean_coverages

    logger.info(
        "Loaded %d BigWigs; sample mean coverage: "
        "min=%.3f, median=%.3f, max=%.3f",
        n_samples,
        np.min(mean_coverages),
        np.median(mean_coverages),
        np.max(mean_coverages),
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

    # Report the selected batch configuration and estimated peak array size.
    estimated_array_mib = (
        window_batch_size * bytes_per_window / 1024**2
    )

    sample_batches = math.ceil(
        n_samples / sample_batch_size
    )

    logger.info(
        "%d BigWigs; %d contigs; %d windows; "
        "%d BigWig batches; "
        "batch limits: %d windows x %d BigWigs; "
        "%d threads; "
        "estimated dense array: %.1f MiB",
        n_samples,
        len(chrom_sizes),
        n_windows,
        sample_batches,
        window_batch_size,
        sample_batch_size,
        args.threads,
        estimated_array_mib,
    )

    output_buffer_size = 8 * 1024 * 1024
    start_time = time.perf_counter()

    bedgraph_path = f"{args.prefix}.normalised-coverage.bedgraph"
    high_mask_path = f"{args.prefix}.high_depth.bed"
    low_mask_path = f"{args.prefix}.low_depth.bed"

    all_cohort_means = np.empty(n_windows, dtype=np.float64)
    processed = 0

    with open(
        bedgraph_path,
        "w",
        buffering=output_buffer_size,
    ) as bedgraph:

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

            batch_end = processed + len(starts)
            all_cohort_means[processed:batch_end] = cohort_mean

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

            processed = batch_end

            elapsed = time.perf_counter() - start_time
            rate = processed / elapsed

            logger.info(
                "Processed %d / %d windows (%.1f%%) at %.0f windows/s",
                processed,
                n_windows,
                100.0 * processed / n_windows,
                rate,
            )
    
    if processed != n_windows:
        raise ValueError(
            f"Processed {processed} windows but expected {n_windows}"
        )

    # Calculate filtering thresholds
    low_threshold, high_threshold = np.percentile(
        all_cohort_means,
        [
            args.low_percentile,
            args.high_percentile,
        ],
    )
    low_selected = all_cohort_means <= low_threshold
    high_selected = all_cohort_means >= high_threshold

    logger.info(
        "Selected %d low-depth windows (%.2f%%) and "
        "%d high-depth windows (%.2f%%)",
        np.count_nonzero(low_selected),
        100.0 * np.mean(low_selected),
        np.count_nonzero(high_selected),
        100.0 * np.mean(high_selected),
    )

    logger.info(
        "Empirical thresholds: %.2f%% = %.8g; %.2f%% = %.8g",
        args.low_percentile,
        low_threshold,
        args.high_percentile,
        high_threshold,
    )

    logger.info(
        "Coverage quantiles: "
        "1%%=%.3f 5%%=%.3f 50%%=%.3f 95%%=%.3f 99%%=%.3f",
        *np.percentile(
            all_cohort_means,
            [1, 5, 50, 95, 99]
        )
    )
    
    write_depth_masks(
        chrom_sizes=chrom_sizes,
        width=args.window_size,
        values=all_cohort_means,
        low_path=low_mask_path,
        high_path=high_mask_path,
        low=low_threshold,
        high=high_threshold,
    )
    del all_cohort_means
    
    elapsed = time.perf_counter() - start_time

    logger.info(
        "Finished: %d windows in %.1f seconds at %.0f windows/s",
        processed,
        elapsed,
        processed / elapsed,
    )


if __name__ == "__main__":
    main()