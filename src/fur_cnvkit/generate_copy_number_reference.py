import argparse
import logging
from pathlib import Path

from fur_cnvkit.normal_vs_normal import perform_normal_vs_normal_comparisons
from fur_cnvkit.utils.cnvkit_utils import (
    run_cnvkit_coverage,
    run_cnvkit_reference,
)
from fur_cnvkit.utils.fur_utils import (
    extract_metadata_files_from_parameter_json,
    split_file_list_by_sample_sex,
)


def configure_logging():
    """
    Define logging configuration
    """
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
    )


def parse_arguments():
    """
    Define and parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate a copy number reference file for a given cohort of samples."
    )
    parser.add_argument(
        "-p",
        "--parameter_file",
        type=Path,
        required=True,
        help="Path to the parameter file.",
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, required=True, help="Path to the output directory."
    )

    return parser.parse_args()


def main():
    # Get command line arguments
    logging.info("Getting command line arguments...")
    args = parse_arguments()

    parameter_file = args.parameter_file
    outdir = args.outdir

    logging.debug(f"Parameter file: {parameter_file}")
    logging.debug(f"Outdir: {outdir}")

    # Extract metadata files from parameter file
    metadata = extract_metadata_files_from_parameter_json(parameter_file)
    normal_bams = metadata["normal_bams"]
    reference_fasta = metadata["reference_fasta"]
    baitset_bed = metadata["baitset_bed"]
    sample_metadata_xlsx = metadata["sample_metadata_xlsx"]
    targets_bed = metadata["targets_bed"]
    antitargets_bed = metadata["antitargets_bed"]

    logging.debug(f"Normal BAMs: {normal_bams}")
    logging.debug(f"Reference FASTA: {str(reference_fasta)}")
    logging.debug(f"Baitset BED: {str(baitset_bed)}")
    logging.debug(f"Sample metadata Excel spreadsheet: {str(sample_metadata_xlsx)}")
    logging.debug(f"Targets BED: {str(targets_bed)}")
    logging.debug(f"Antitargets BED: {str(antitargets_bed)}")

    # Process normal BAMs to get a dictionary of BAMs separated by sex
    logging.info("Separating normal BAMs based on their sex ...")
    sex_separated_normal_bams = split_file_list_by_sample_sex(
        normal_bams, sample_metadata_xlsx
    )
    logging.debug(f"Split normal BAM files by sex: {sex_separated_normal_bams}")

    # Use CNVKit to generate a copy number reference for each sex
    for sex in sex_separated_normal_bams:
        logging.info(
            f"Running CNVKit copy number reference generation for {sex} samples ..."
        )

        sex_outdir = outdir / sex
        sex_outdir.mkdir(parents=True, exist_ok=True)

        logging.debug(
            f"Output directory for {sex} samples generated at {str(sex_outdir)}"
        )

        # Extract the sex-separated normal BAMs for the given sex from the dictionary
        sex_normal_bams = sex_separated_normal_bams[sex]

        # Initialise a list to store coverage files for this sex
        sex_normal_coverage_files = []

        # Make a directory for coverage files, if it does not exist
        coverage_file_dir = sex_outdir / "coverage_files"
        coverage_file_dir.mkdir(parents=True, exist_ok=True)

        # Loop through the normal sample BAMs for this sex
        for sample_bam in sex_normal_bams:
            # Run cnvkit.py coverage for this sample on target regions
            sample_target_coverage_file = run_cnvkit_coverage(
                sample_bam, targets_bed, coverage_file_dir
            )

            # Run cnvkit.py coverage for this sample on antitarget regions
            sample_antitarget_coverage_file = run_cnvkit_coverage(
                sample_bam, antitargets_bed, coverage_file_dir
            )

            # Append the coverage files for this sample to the list
            sex_normal_coverage_files.extend(
                [sample_target_coverage_file, sample_antitarget_coverage_file]
            )

        logging.debug(f"Coverage files for {sex} samples: {sex_normal_coverage_files}")

        # Perform pairwise normal vs. normal comparisons to filter coverage files
        logging.info(f"Performing normal vs. normal comparisons for {sex} samples ...")

        # Make a directory for the normal vs. normal comparisons, if it does not exist yet
        normal_vs_normal_dir = sex_outdir / "normal_vs_normal"
        normal_vs_normal_dir.mkdir(parents=True, exist_ok=True)

        logging.debug(
            f"Normal vs. normal comparison directory for {sex} samples generated at {str(normal_vs_normal_dir)}"
        )

        # Perform normal vs. normal comparisons to produce a list of filtered coverage files
        logging.info(
            "Performing normal vs. normal comparisons to filter out potentially spurious coverage files ..."
        )

        filtered_coverage_files = perform_normal_vs_normal_comparisons(
            sex_normal_coverage_files,
            reference_fasta,
            sample_metadata_xlsx,
            normal_vs_normal_dir,
        )

        logging.debug(f"Filtered coverage files: {filtered_coverage_files}")

        # Run CNVKit reference to create a new copy number reference using the filtered coverage files
        logging.info(
            f"Generating a copy number reference for {sex} samples using the filtered coverage files ..."
        )

        sex_copy_number_reference_file = run_cnvkit_reference(
            coverage_files=filtered_coverage_files,
            reference_fasta=reference_fasta,
            output_prefix=sex,
            outdir=sex_outdir,
            sex=sex,
        )

        logging.debug(
            f"Copy number reference for {sex} generated at {str(sex_copy_number_reference_file)}"
        )


if __name__ == "__main__":
    configure_logging()
    main()
