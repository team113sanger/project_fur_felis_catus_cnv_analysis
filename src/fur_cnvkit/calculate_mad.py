import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Optional
import math

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from fur_cnvkit.utils.fur_utils import get_sample_id_from_file_path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run MAD calculation pipeline on CNVkit files."
    )
    parser.add_argument(
        "--cnr_files",
        nargs="+",
        type=Path,
        required=True,
        help="List of CNVkit .cnr files.",
    )
    parser.add_argument(
        "--cns_files",
        nargs="+",
        type=Path,
        help="List of corresponding CNVkit .cns files.",
        required=False,
    )
    parser.add_argument(
        "--outdir", type=Path, required=True, help="Output directory for results."
    )
    parser.add_argument(
        "--prefix", type=str, default="result", help="Prefix for output files."
    )
    parser.add_argument(
        "--log2_column",
        type=str,
        default="log2",
        help="Name of the log2 column (default 'log2').",
    )
    parser.add_argument(
        "--show_plot", action="store_true", help="Display the plot interactively."
    )
    parser.add_argument(
        "--z_threshold",
        type=float,
        default=3.5,
        help="Modified z-score threshold for filtering noisy samples (default 3.5).",
    )
    return parser.parse_args()


def calculate_mad(cnr_file: Path, log2_column: str = "log2") -> float:
    """
    Calculate the Median Absolute Deviation (MAD) for a given CNVkit .cnr file.
    """
    try:
        data = pd.read_csv(cnr_file, sep="\t")
    except Exception as e:
        logging.error(f"Error reading file {cnr_file}: {e}")
        raise

    if log2_column not in data.columns:
        raise ValueError(f"Column '{log2_column}' not found in file {cnr_file}")

    log2_values = data[log2_column]
    median_value = log2_values.median()
    mad_value = (log2_values - median_value).abs().median()
    return mad_value


def calculate_average_segment_size(cns_file: Path) -> Optional[float]:
    """
    Calculate the average segment size (in base pairs) from a CNVkit .cns file.
    If the file does not exist or is invalid, returns None.
    """
    if not cns_file.exists():
        logging.warning(f"Segment file not found: {cns_file}")
        return None

    try:
        seg_data = pd.read_csv(cns_file, sep="\t")
    except Exception as e:
        logging.error(f"Error reading file {cns_file}: {e}")
        return None

    required_cols = {"chromosome", "start", "end"}
    if not required_cols.issubset(seg_data.columns):
        logging.warning(
            f"Columns {required_cols} not found in {cns_file}; cannot compute segment size."
        )
        return None

    seg_data["length"] = seg_data["end"] - seg_data["start"]
    avg_size = seg_data["length"].mean()
    return avg_size


def adjusted_metric(mad_value: float, avg_segment_size: Optional[float]) -> float:
    """
    Compute an adjusted MAD metric that combines MAD and average segment size.
    Large segments reduce the effective score.

    Formula:
        AdjustedMetric = MAD / log2(avg_segment_size + 2)
    If avg_segment_size is None or invalid, simply return MAD.
    """
    if avg_segment_size is None or avg_segment_size <= 0:
        return mad_value
    return mad_value / math.log2(avg_segment_size + 2)


def filter_noisy_samples_by_zscore(
    results_df: pd.DataFrame, metric_col: str, z_threshold: float = 3.5
) -> List[str]:
    """
    Filter out samples based on a modified z-score of the chosen metric.
    The modified z-score is calculated as:

        modified_z = 0.6745 * (x - median) / MAD_metric

    Samples with a modified z-score greater than z_threshold are considered noisy.
    """
    metric_values = results_df[metric_col]
    median_val = metric_values.median()
    mad_metric = (metric_values - median_val).abs().median()
    if mad_metric == 0:
        logging.warning("MAD of metric is zero; no outlier filtering applied.")
        return []

    # Compute the modified z-score for each sample.
    modified_z = 0.6745 * (metric_values - median_val) / mad_metric
    # Filter samples where the modified z-score exceeds the threshold.
    noisy_samples = results_df.loc[modified_z > z_threshold, "Sample"].tolist()
    logging.info(
        f"Filtering using z_threshold {z_threshold}. Median: {median_val:.4f}, MAD: {mad_metric:.4f}"
    )
    logging.info(f"Noisy samples (modified z-score > {z_threshold}): {noisy_samples}")
    return noisy_samples


def plot_metric_distribution(
    results_df: pd.DataFrame,
    metric_col: str,
    outdir: Path,
    prefix: str,
    show_plot: bool = False,
) -> None:
    """
    Plot and save the distribution of a chosen metric (e.g. MAD or AdjustedMADMetric).
    """
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(data=results_df, x=metric_col, bins=20, kde=True)
    plt.title(f"Distribution of {metric_col} across Samples")
    plt.xlabel(metric_col)
    plt.ylabel("Frequency")

    outdir.mkdir(parents=True, exist_ok=True)
    plot_path = outdir / f"{prefix}.{metric_col}_distribution.png"
    plt.savefig(plot_path)
    if show_plot:
        plt.show()
    else:
        plt.close()
    logging.info(f"{metric_col} distribution plot saved to {plot_path}")


def build_cns_mapping(cns_files: Optional[List[Path]]) -> dict:
    """
    Build a mapping from sample ID to corresponding .cns file.
    """
    mapping = {}
    if cns_files:
        for cns in cns_files:
            sample_id = get_sample_id_from_file_path(cns)
            mapping[sample_id] = cns
    return mapping


def process_sample(
    cnr_file: Path, cns_mapping: dict, log2_column: str
) -> Optional[dict]:
    """
    Process a single sample: calculate MAD, determine the corresponding .cns file,
    compute average segment size, and calculate the adjusted MAD metric.
    """
    sample_id = get_sample_id_from_file_path(cnr_file)
    try:
        mad_value = calculate_mad(cnr_file, log2_column=log2_column)
    except Exception as e:
        logging.error(
            f"Error calculating MAD for sample {sample_id} from file {cnr_file}: {e}"
        )
        return None

    # Determine the corresponding .cns file.
    cns_file = cns_mapping.get(sample_id, cnr_file.with_suffix(".cns"))
    avg_seg_size = calculate_average_segment_size(cns_file)
    adj_val = adjusted_metric(mad_value, avg_seg_size)

    return {
        "Sample": sample_id,
        "MAD": mad_value,
        "AvgSegSize": avg_seg_size,
        "AdjustedMADMetric": adj_val,
        "CNR": str(cnr_file),
        "CNS": str(cns_file) if cns_file.exists() else None,
    }


def save_results_df(results_df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """
    Save the DataFrame of results to a TSV file.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / f"{prefix}.metrics.tsv"
    results_df.to_csv(tsv_path, sep="\t", index=False)
    logging.info(f"Metrics (including MAD and adjusted MAD metric) saved to {tsv_path}")


def generate_plots(
    results_df: pd.DataFrame, outdir: Path, prefix: str, show_plot: bool
) -> None:
    """
    Generate and save distribution plots for the chosen metrics.
    """
    plot_metric_distribution(
        results_df, metric_col="MAD", outdir=outdir, prefix=prefix, show_plot=show_plot
    )
    plot_metric_distribution(
        results_df,
        metric_col="AdjustedMADMetric",
        outdir=outdir,
        prefix=prefix,
        show_plot=show_plot,
    )


def run_mad_calculation_pipeline(
    cnr_files: List[Path],
    outdir: Path,
    prefix: str,
    log2_column: str = "log2",
    show_plot: bool = False,
    cns_files: Optional[List[Path]] = None,
    z_threshold: float = 3.5,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Run the pipeline by processing each CNVkit file:
      1. Build mapping from sample ID to .cns file (if provided).
      2. Process each .cnr file.
      3. Save the results and generate plots.
      4. Filter out noisy samples based on the modified z-score.
    """
    # Build mapping for .cns files.
    cns_mapping = build_cns_mapping(cns_files)

    # Process each sample.
    results = []
    for cnr_file in cnr_files:
        sample_result = process_sample(cnr_file, cns_mapping, log2_column)
        if sample_result is not None:
            results.append(sample_result)

    results_df = pd.DataFrame(results)

    # Save results and generate plots.
    save_results_df(results_df, outdir, prefix)
    generate_plots(results_df, outdir, prefix, show_plot)

    # Filter out noisy samples using modified z-score.
    filtered_sample_ids = filter_noisy_samples_by_zscore(
        results_df, metric_col="AdjustedMADMetric", z_threshold=z_threshold
    )
    return filtered_sample_ids, results_df


def main():
    args = parse_args()

    filtered_samples, results_df = run_mad_calculation_pipeline(
        cnr_files=args.cnr_files,
        outdir=args.outdir,
        prefix=args.prefix,
        log2_column=args.log2_column,
        show_plot=args.show_plot,
        cns_files=args.cns_files,
        z_threshold=args.z_threshold,
    )
    # Print filtered sample IDs to stdout.
    print(f"Samples identified as noisy (modified z-score > {args.z_threshold}):")
    for sample in filtered_samples:
        print(sample)


if __name__ == "__main__":
    main()
