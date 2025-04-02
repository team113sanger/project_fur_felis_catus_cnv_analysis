import logging
from pathlib import Path
import typing as t

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

from fur_cnvkit.utils.cnvkit_utils import (
    validate_sample_coverage_files,
    run_cnvkit_reference,
    run_cnvkit_fix,
    run_cnvkit_genemetrics,
    parse_genemetrics_file,
)
from fur_cnvkit.utils.fur_utils import (
    filter_files_by_exclude_samples,
    get_sample_ids_for_file_list,
    get_sample_specific_files,
    map_sample_ids_to_study_ids,
    get_sample_sex,
)

# Set up logging
logger = logging.getLogger(__name__)


def perform_normal_vs_normal_comparisons(
    normal_coverage_files: t.List[Path],
    reference_fasta: Path,
    sample_metadata_xlsx: Path,
    outdir: Path,
):
    logger.info("Performing normal vs normal comparisons ...")
    logger.debug(f"Normal coverage files: {normal_coverage_files}")

    # Get normal sample IDs
    logger.info("Getting normal sample IDs for normal coverage files ...")
    normal_sample_ids = get_sample_ids_for_file_list(normal_coverage_files)

    logger.debug(f"Normal sample IDs: {normal_sample_ids}")

    # Map normal sample IDs to their corresponding study ID
    logger.info("Mapping normal sample IDs to study IDs ...")
    sample_study_dict = map_sample_ids_to_study_ids(
        sample_ids=normal_sample_ids, sample_metadata_xlsx=sample_metadata_xlsx
    )

    logger.debug(f"Sample study dictionary: {sample_study_dict}")

    # Create a dictionary of collated log2 ratios CSV files for each study
    # Key: study ID
    # Value: collated log2 ratios CSV file
    study_collated_log2_ratios_csv_dict = {}

    # Loop through each study and perform normal vs normal comparisons
    for study_id, sample_ids in sample_study_dict.items():
        logger.info(f"Performing normal vs normal comparisons for study {study_id} ...")

        # Create study-specific output directory
        study_outdir = outdir / study_id
        study_outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Study-specific output directory: {study_outdir}")

        # Get sample-specific normal coverage files for the study
        study_normal_coverage_files_dict = {}

        for sample_id in sample_ids:
            study_normal_coverage_files_dict[sample_id] = get_sample_specific_files(
                files=normal_coverage_files,
                sample=sample_id,
                suffixes=[".targetcoverage.cnn", ".antitargetcoverage.cnn"],
            )

        logger.debug(
            f"Study normal coverage files dictionary: {study_normal_coverage_files_dict}"
        )

        # Perform normal vs normal comparisons
        study_collated_log2_ratios_csv = perform_normal_vs_normal_comparisons_for_study(
            study_id=study_id,
            sample_ids=sample_ids,
            study_normal_coverage_files_dict=study_normal_coverage_files_dict,
            reference_fasta=reference_fasta,
            sample_metadata_xlsx=sample_metadata_xlsx,
            outdir=study_outdir,
        )

        # Add the collated log2 ratios CSV file to the dictionary
        study_collated_log2_ratios_csv_dict[study_id] = study_collated_log2_ratios_csv

    # After all studies have been processed, collate the average median log2(FC) values across studies
    combined_median_of_median_log2_ratios_df = (
        collate_median_of_median_log2_ratios_across_studies(
            study_collated_log2_ratios_csv_dict=study_collated_log2_ratios_csv_dict,
            outdir=outdir,
        )
    )

    # Filter the combined DataFrame to only the top and bottom 10% of samples per study
    top_bottom_10pct_df = get_top_bottom_10pct_samples(
        combined_median_of_median_log2_ratios_df
    )

    # Get a list of the top and bottom 10% sample IDs
    top_10pct_samples = top_bottom_10pct_df[top_bottom_10pct_df["Group"] == "top_10pct"]
    bottom_10pct_samples = top_bottom_10pct_df[
        top_bottom_10pct_df["Group"] == "bottom_10pct"
    ]
    sample_ids_10pct = list(top_10pct_samples["Reference Sample"]) + list(
        bottom_10pct_samples["Reference Sample"]
    )

    logger.debug(
        f"Samples filtered out by normal vs. normal comparisons: {sample_ids_10pct}"
    )

    # Filter out the top and bottom 10% samples from the normal coverage files
    filtered_normal_coverage_files = filter_files_by_exclude_samples(
        normal_coverage_files, sample_ids_10pct
    )

    return filtered_normal_coverage_files


def get_top_bottom_10pct_samples(
    combined_df: pd.DataFrame, lower_quantile: float = 0.1, upper_quantile: float = 0.9
) -> pd.DataFrame:
    """
    Given a DataFrame with columns ['Study ID', 'Reference Sample', 'Average Median Log2(FC)'],
    return a DataFrame of only the samples in the bottom and top 10% (by default) per study.

    :param combined_df: DataFrame containing the merged results across studies
    :param lower_quantile: The lower quantile cutoff, default 0.10 (10%)
    :param upper_quantile: The upper quantile cutoff, default 0.90 (90%)
    :return: A DataFrame filtered to only the top/bottom 10% of samples for each study,
             with a new column 'Group' labeling them as 'bottom_10pct' or 'top_10pct'.
    """

    # Ensure the required columns exist
    required_columns = {"Study ID", "Reference Sample", "Average Median Log2(FC)"}
    missing_cols = required_columns - set(combined_df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # List to accumulate bottom/top results per study
    filtered_dfs = []

    # Group by Study ID
    for study_id, group_df in combined_df.groupby("Study ID"):
        # Calculate the cutoff values for this study
        low_cutoff = group_df["Average Median Log2(FC)"].quantile(lower_quantile)
        high_cutoff = group_df["Average Median Log2(FC)"].quantile(upper_quantile)

        # Bottom 10% samples (those below the lower_quantile)
        bottom_10pct_df = group_df[
            group_df["Average Median Log2(FC)"] < low_cutoff
        ].copy()
        bottom_10pct_df["Group"] = "bottom_10pct"

        # Top 10% samples (those above the upper_quantile)
        top_10pct_df = group_df[
            group_df["Average Median Log2(FC)"] > high_cutoff
        ].copy()
        top_10pct_df["Group"] = "top_10pct"

        # Append both subsets to the master list
        filtered_dfs.extend([bottom_10pct_df, top_10pct_df])

    # Concatenate all filtered DataFrames
    if filtered_dfs:
        result_df = pd.concat(filtered_dfs, ignore_index=True)
    else:
        # If no data matched the criteria (corner case), return empty
        result_df = pd.DataFrame(
            columns=["Study ID", "Reference Sample", "Average Median Log2(FC)", "Group"]
        )

    return result_df


def collate_median_of_median_log2_ratios_across_studies(
    study_collated_log2_ratios_csv_dict: t.Dict[str, Path],
    outdir: Path,
    output_filename: str = "combined_average_median_log2_fc.csv",
) -> pd.DataFrame:
    """
    Given a dictionary mapping study IDs to the path of their 'collated_log2_ratios.csv',
    compute the average median log2(FC) per reference sample in each study and combine
    into a single DataFrame with an added 'Study ID' column.

    :param study_collated_log2_ratios_csv_dict: Dictionary of {study_id: Path_to_collated_log2_ratios.csv}
    :param outdir: Output directory where the combined CSV will be written
    :param output_filename: Name of the output CSV file
    :return: A combined pandas DataFrame containing columns:
             ['Study ID', 'Reference Sample', 'Average Median Log2(FC)']
    """
    combined_df_list = []

    for study_id, log2_csv_path in study_collated_log2_ratios_csv_dict.items():
        logger.info(f"Processing study {study_id} using file {log2_csv_path}")

        # Load the collated log2 data
        log2_data = pd.read_csv(log2_csv_path)

        # Calculate the median Log2(FC) for each (Reference Sample, Comparator Sample)
        median_per_comparator = (
            log2_data.groupby(["Reference Sample", "Comparator Sample"])["Log2(FC)"]
            .median()
            .reset_index()
        )

        average_median_log2 = (
            median_per_comparator.groupby("Reference Sample")["Log2(FC)"]
            .median()
            .reset_index()
            .rename(columns={"Log2(FC)": "Average Median Log2(FC)"})
        )

        # Add the 'Study ID' column
        average_median_log2["Study ID"] = study_id

        # Collect this into our list of DataFrames
        combined_df_list.append(average_median_log2)

    # Concatenate all study-specific DataFrames
    if not combined_df_list:
        logger.warning("No data found. Returning an empty DataFrame.")
        return pd.DataFrame()

    combined_df = pd.concat(combined_df_list, ignore_index=True)
    combined_df = combined_df[
        ["Study ID", "Reference Sample", "Average Median Log2(FC)"]
    ]

    # Write out the combined CSV
    output_csv_path = outdir / output_filename
    combined_df.to_csv(output_csv_path, index=False)
    logger.info(f"Combined average median log2(FC) CSV written to: {output_csv_path}")

    # Create a boxplot of the average median log2(FC) values
    create_average_median_log2_fc_boxplot(combined_df, outdir)

    return combined_df


def create_average_median_log2_fc_boxplot(
    average_median_log2_df: pd.DataFrame, outdir: Path
):
    # Create the plot
    sns.boxplot(
        data=average_median_log2_df,
        x="Average Median Log2(FC)",
        y="Study ID",
        hue="Study ID",
    )

    # Save the plot as a PDF file
    output_plot = outdir / "average_median_log2_fc_boxplot.pdf"
    plt.savefig(output_plot)


def calculate_median_of_median_log2_ratios(input_log2_csv: Path, output_csv: Path):
    """Calculate the average median log2(FC) value for each reference sample in a pairwise normal vs. normal comparison"""

    # Load the data
    log2_data = pd.read_csv(input_log2_csv)

    # Group by 'Reference Sample' and 'Comparator Sample' to calculate the median Log2(FC) for each comparator sample
    median_per_comparator = (
        log2_data.groupby(["Reference Sample", "Comparator Sample"])["Log2(FC)"]
        .median()
        .reset_index()
    )

    # Group by 'Reference Sample' again to calculate the average of these medians for each reference sample
    average_median_log2 = (
        median_per_comparator.groupby("Reference Sample")["Log2(FC)"]
        .median()
        .reset_index()
    )

    # Rename the columns for clarity
    average_median_log2.columns = ["Reference Sample", "Average Median Log2(FC)"]

    average_median_log2.to_csv(output_csv)


def perform_normal_vs_normal_comparisons_for_study(
    study_id: str,
    sample_ids: t.List[str],
    study_normal_coverage_files_dict: t.Dict[str, t.List[Path]],
    reference_fasta: Path,
    sample_metadata_xlsx: Path,
    outdir: Path,
):
    # Create a dictionary of genemetrics files for each reference sample
    # Key: reference sample ID
    # Value: list of genemetrics files for comparison samples when using the reference sample as the copy number reference
    reference_vs_comparison_genemetrics_files_dict = {}

    # Loop through each sample and perform normal vs normal comparisons
    for reference_sample_id in sample_ids:
        logger.info(
            f"Performing normal vs normal comparisons using sample {reference_sample_id} as the reference sample..."
        )

        # Create sample-specific output directory
        reference_sample_outdir = outdir / reference_sample_id
        reference_sample_outdir.mkdir(parents=True, exist_ok=True)

        # Perform normal vs normal comparisons for the sample
        comparison_sample_genemetrics_files = (
            perform_normal_vs_normal_comparisons_for_sample(
                reference_sample_id=reference_sample_id,
                study_normal_coverage_files_dict=study_normal_coverage_files_dict,
                reference_fasta=reference_fasta,
                sample_metadata_xlsx=sample_metadata_xlsx,
                outdir=reference_sample_outdir,
            )
        )

        # Add the genemetrics files to the dictionary
        reference_vs_comparison_genemetrics_files_dict[
            reference_sample_id
        ] = comparison_sample_genemetrics_files

    # Collate the log2 ratios for each gene for each comparison sample
    collated_log2_ratios_csv = outdir / "collated_log2_ratios.csv"
    collated_log2_ratios_df = collate_log2_ratios(
        reference_vs_comparison_genemetrics_files_dict, collated_log2_ratios_csv
    )

    # Create a boxplot of the log2 ratios for each gene from each comparison
    create_log2_ratio_boxplot(collated_log2_ratios_df, outdir)

    return collated_log2_ratios_csv


def collate_log2_ratios(
    reference_vs_comparison_genemetrics_files_dict: t.Dict[str, t.List[Path]],
    output_csv: Path,
):
    logger.info("Collating log2 ratios for each gene for each comparison ...")

    log2_data = []

    for (
        genemetircs_file_list
    ) in reference_vs_comparison_genemetrics_files_dict.values():
        for genemetrics_file in genemetircs_file_list:
            comparison_sample_id = genemetrics_file.stem.split("_vs_")[0]
            reference_sample_id = genemetrics_file.stem.split("_vs_")[1].split(".")[0]

            log2_values = parse_genemetrics_file(genemetrics_file)

            for value in log2_values:
                log2_data.append(
                    {
                        "Reference Sample": reference_sample_id,
                        "Comparator Sample": comparison_sample_id,
                        "Log2(FC)": value,
                    }
                )

    # Convert the log2 data to a DataFrame
    log2_df = pd.DataFrame(log2_data)

    # Write the log2 data to a CSV file
    log2_df.to_csv(output_csv, index=False)
    logger.info(f"Log2 ratios collated and written to {output_csv}")

    return log2_df


def create_log2_ratio_boxplot(log2_ratios_df: pd.DataFrame, outdir: Path):
    logger.info(
        "Creating a boxplot of the log2 ratios for each gene for each comparison ..."
    )

    # Combine 'Reference Sample' and 'Comparator Sample' into a single column
    log2_ratios_df["Combined Sample"] = (
        log2_ratios_df["Reference Sample"]
        + " vs "
        + log2_ratios_df["Comparator Sample"]
    )

    # Sort the DataFrame by 'Combined Sample' and 'Reference Sample'
    log2_ratios_df = log2_ratios_df.sort_values(
        by=["Combined Sample", "Reference Sample"]
    )

    # Get the unique reference samples and sort them alphabetically
    unique_references = sorted(log2_ratios_df["Reference Sample"].unique())
    num_reference_samples = len(unique_references)

    # Create subplots for each reference sample
    fig, axs = plt.subplots(
        num_reference_samples, 1, figsize=(15, num_reference_samples * 5), squeeze=False
    )

    # Generate subplots for each reference sample
    for i, ref_sample in enumerate(unique_references):
        ax = axs[i, 0]
        subset = log2_ratios_df[log2_ratios_df["Reference Sample"] == ref_sample]

        sns.boxplot(x="Combined Sample", y="Log2(FC)", data=subset, ax=ax)
        ax.set_title(f"Log2(FC) Distribution for {ref_sample}")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_xlabel("Combined Sample")
        ax.set_ylabel("Log2(FC)")
        ax.legend().set_visible(False)

        # Add horizontal dotted lines
        for y_value in range(-1, 2):
            ax.axhline(y=y_value, color="grey", linestyle="--", linewidth=0.5)

        # Mirror y-axis labels on the right side
        ax_right = ax.twinx()
        ax_right.set_ylabel("Log2(FC)")
        ax_right.set_yticks(ax.get_yticks())
        ax_right.set_ylim(ax.get_ylim())
        ax_right.set_yticklabels(ax.get_yticklabels())

    plt.tight_layout()

    # Save the plot as a PDF file
    output_plot = outdir / "log2_ratio_boxplot.pdf"
    with PdfPages(output_plot) as pdf:
        pdf.savefig(fig)
        plt.close(fig)

    logging.info(f"Boxplot saved to {str(output_plot)}")


def perform_normal_vs_normal_comparisons_for_sample(
    reference_sample_id: str,
    study_normal_coverage_files_dict: t.Dict[str, t.List[Path]],
    reference_fasta: Path,
    sample_metadata_xlsx: Path,
    outdir: Path,
):
    logger.info(
        f"Performing normal vs normal comparisons for sample {reference_sample_id} ..."
    )

    # Perform normal vs normal comparisons for the reference sample

    # 1. Get the sex of the reference sample
    logger.info("Getting the sex for the reference sample ...")
    reference_sample_sex = get_sample_sex(
        reference_sample_id, sample_metadata_xlsx=sample_metadata_xlsx
    )

    logger.debug(f"Sex for {reference_sample_id}: {reference_sample_sex}")

    # 2. Check that the given sample has all expected coverage files
    logger.info(
        "Checking that the reference sample has all expected coverage files ..."
    )
    (
        reference_sample_target_coverage_file,
        reference_sample_antitarget_coverage_file,
    ) = validate_sample_coverage_files(
        reference_sample_id, study_normal_coverage_files_dict[reference_sample_id]
    )

    # 3. Create a copy number reference for the sample
    logger.info("Creating a copy number reference for the sample ...")
    reference_sample_copy_number_reference_file = run_cnvkit_reference(
        coverage_files=[
            reference_sample_target_coverage_file,
            reference_sample_antitarget_coverage_file,
        ],
        reference_fasta=reference_fasta,
        output_prefix=reference_sample_id,
        outdir=outdir,
        sex=reference_sample_sex,
    )

    logger.info(
        f"Copy number reference file created for sample {reference_sample_id}: {reference_sample_copy_number_reference_file}"
    )

    # 4. Perform normal vs normal comparisons for the sample
    comparison_sample_genetrics_files = compare_all_other_study_samples_to_reference_sample(
        reference_sample_id=reference_sample_id,
        reference_sample_copy_number_reference_file=reference_sample_copy_number_reference_file,
        study_normal_coverage_files_dict=study_normal_coverage_files_dict,
        sample_metadata_xlsx=sample_metadata_xlsx,
        outdir=outdir,
    )

    # 5. Return the genemetrics files for the comparison samples
    return comparison_sample_genetrics_files


def compare_all_other_study_samples_to_reference_sample(
    reference_sample_id: str,
    reference_sample_copy_number_reference_file: Path,
    study_normal_coverage_files_dict: t.Dict[str, t.List[Path]],
    sample_metadata_xlsx: Path,
    outdir: Path,
):
    # Get a set of every other normal sample ID in the study to compare to the reference sample
    def get_other_sample_ids(sample_ids):
        return set(sample_ids) - {reference_sample_id}

    other_sample_ids = get_other_sample_ids(study_normal_coverage_files_dict.keys())

    # Create a list of genemetircs files for each comparison sample
    genemetrics_files = []

    # Loop through each other sample and compare to the reference sample
    for comparison_sample_id in other_sample_ids:
        logger.info(
            f"Comparing sample {comparison_sample_id} to reference sample {reference_sample_id} ..."
        )

        # Get the coverage files for the comparison sample
        comparison_sample_coverage_files = study_normal_coverage_files_dict[
            comparison_sample_id
        ]
        (
            comparison_sample_target_coverage_file,
            comparison_sample_antitarget_coverage_file,
        ) = validate_sample_coverage_files(
            comparison_sample_id, comparison_sample_coverage_files
        )

        # Run cnvkit.py fix to adjust other sample's (anti)target coverage files to the reference sample's copy number reference
        comparison_sample_ratio_file = run_cnvkit_fix(
            comparison_sample_target_coverage_file,
            comparison_sample_antitarget_coverage_file,
            reference_sample_copy_number_reference_file,
            comparison_sample_id,
            outdir,
        )

        # Get the sex of the comparison sample
        comparison_sample_sex = get_sample_sex(
            comparison_sample_id, sample_metadata_xlsx=sample_metadata_xlsx
        )

        # Run CNVKit genemetrics to compare the comparison sample to the reference sample
        comparison_id = f"{comparison_sample_id}_vs_{reference_sample_id}"
        comparison_sample_genemetrics_file = run_cnvkit_genemetrics(
            ratio_file=comparison_sample_ratio_file,
            threshold=0.00001,  # Set low so all genes are reported in output
            min_probes=3,  # Set low so all genes are reported in output
            output_prefix=comparison_id,
            outdir=outdir,
            sex=comparison_sample_sex,
        )

        # Add the genemetircs file to the list
        genemetrics_files.append(comparison_sample_genemetrics_file)

    return genemetrics_files
