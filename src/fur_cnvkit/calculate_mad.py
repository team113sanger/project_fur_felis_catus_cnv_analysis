import argparse
import logging
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from fur_cnvkit.utils.fur_utils import get_sample_id_from_file_path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def calculate_mad(cnr_file: Path, log2_column: str = "log2") -> float:
    """
    Calculate the Median Absolute Deviation (MAD) for a given CNVkit .cnr file.

    Parameters:
      cnr_file (Path): Path to the CNVkit .cnr file.
      log2_column (str): Column name containing log2 fold-change values (default "log2").

    Returns:
      float: The calculated MAD value.
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


def filter_top_mad_samples(
    results_df: pd.DataFrame, top_percent: float = 0.2
) -> List[str]:
    """
    Filter out samples with the top `top_percent` highest MAD values.

    Parameters:
      results_df (pd.DataFrame): DataFrame with at least 'Sample' and 'MAD' columns.
      top_percent (float): Proportion of samples to filter out (default is 0.2, i.e. top 20%).

    Returns:
      List[str]: List of sample IDs that are in the top `top_percent` highest MAD values.
    """
    threshold = results_df["MAD"].quantile(1 - top_percent)
    logging.info(f"MAD threshold for top {top_percent*100:.0f}%: {threshold}")
    filtered_samples = results_df[results_df["MAD"] >= threshold]["Sample"].tolist()
    return filtered_samples


def plot_mad_distribution(
    results_df: pd.DataFrame, outdir: Path, prefix: str, show_plot: bool = False
) -> None:
    """
    Plot and save the distribution of MAD values.

    Parameters:
      results_df (pd.DataFrame): DataFrame containing the MAD results.
      outdir (Path): Output directory.
      prefix (str): Prefix for the saved plot file.
      show_plot (bool): If True, display the plot interactively.
    """
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(data=results_df, x="MAD", bins=20, kde=True)
    plt.title("Distribution of MAD across Samples")
    plt.xlabel("Median Absolute Deviation (MAD)")
    plt.ylabel("Frequency")

    outdir.mkdir(parents=True, exist_ok=True)
    plot_path = outdir / f"{prefix}.mad_distribution.png"
    plt.savefig(plot_path)
    if show_plot:
        plt.show()
    else:
        plt.close()
    logging.info(f"MAD distribution plot saved to {plot_path}")


def run_mad_calculation_pipeline(
    cnr_files: List[Path],
    outdir: Path,
    prefix: str,
    top_percent: float = 0.2,
    log2_column: str = "log2",
) -> Tuple[List[str], pd.DataFrame]:
    """
    Run the complete MAD calculation pipeline:
      - Calculate MAD for each provided .cnr file.
      - Save the MAD results to a TSV file.
      - Generate a MAD distribution plot.
      - Filter out the samples with the top `top_percent` highest MAD values.

    Parameters:
      cnr_files (List[Path]): List of paths to the .cnr files.
      outdir (Path): Output directory for results.
      prefix (str): Prefix for output files.
      top_percent (float): Proportion of samples to filter out (default is 0.2 for top 20%).
      log2_column (str): Name of the log2 column (default "log2").

    Returns:
      Tuple[List[str], pd.DataFrame]:
        - List of sample IDs (as strings) that are in the top `top_percent` (to be filtered out).
        - DataFrame containing the MAD results.
    """
    results = []
    for cnr_file in cnr_files:
        sample_id = get_sample_id_from_file_path(cnr_file)
        try:
            mad_value = calculate_mad(cnr_file, log2_column=log2_column)
        except Exception as e:
            logging.error(
                f"Error calculating MAD for sample {sample_id} from file {cnr_file}: {e}"
            )
            continue
        results.append({"Sample": sample_id, "MAD": mad_value, "CNR": str(cnr_file)})
    results_df = pd.DataFrame(results)

    # Save the MAD results to a TSV file.
    tsv_path = outdir / f"{prefix}.mad.tsv"
    results_df.to_csv(tsv_path, sep="\t", index=False)
    logging.info(f"MAD results saved to {tsv_path}")

    # Generate a MAD distribution plot.
    plot_mad_distribution(results_df, outdir, prefix=prefix, show_plot=False)

    # Filter out samples with the top MAD values.
    filtered_sample_ids = filter_top_mad_samples(results_df, top_percent=top_percent)
    logging.info(
        f"Samples to filter out (top {top_percent*100:.0f}% highest MAD): {filtered_sample_ids}"
    )

    return filtered_sample_ids, results_df


def main():
    parser = argparse.ArgumentParser(
        description="Run MAD calculation on CNVkit .cnr files."
    )
    parser.add_argument(
        "cnr_files",
        nargs="+",
        type=Path,
        help="List of CNVkit .cnr files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for MAD results.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="result",
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--top_percent",
        type=float,
        default=0.2,
        help="Proportion of samples to filter out (default 0.2 for top 20%).",
    )
    parser.add_argument(
        "--log2_column",
        type=str,
        default="log2",
        help="Name of the log2 column (default 'log2').",
    )
    parser.add_argument(
        "--show_plot",
        action="store_true",
        help="Display the plot interactively.",
    )
    args = parser.parse_args()

    filtered_samples, results_df = run_mad_calculation_pipeline(
        cnr_files=args.cnr_files,
        outdir=args.outdir,
        prefix=args.prefix,
        top_percent=args.top_percent,
        log2_column=args.log2_column,
    )
    # Print filtered sample IDs to stdout.
    print(f"Samples filtered out (top {args.top_percent*100:.0f}% highest MAD):")
    for sample in filtered_samples:
        print(sample)


if __name__ == "__main__":
    main()
