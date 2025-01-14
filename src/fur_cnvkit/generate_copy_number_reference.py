import argparse
import logging
from pathlib import Path

# import typing as t

from utils.cnvkit_utils import run_cnvkit_coverage

# from utils.file_format_checker import validate_bam_files
from utils.fur_utils import (
    extract_metadata_files_from_config_json,
    # remove_unwanted_sample_files,
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
    # parser.add_argument(
    #     "-n",
    #     "--normal_bams",
    #     type=Path,
    #     nargs="+",
    #     required=True,
    #     help="Path to normal BAM(s) to generate copy number referemce with.",
    # )
    # parser.add_argument(
    #     "-e",
    #     "--exclude_file",
    #     type=Path,
    #     required=False,
    #     help="Exclude file containing samples to exclude from copy number reference (One sample on each line).",
    # )
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        required=True,
        help="Path to the configuration file.",
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, required=True, help="Path to the output directory."
    )

    return parser.parse_args()


# def process_normal_bams(
#     normal_bams: t.List[Path],
#     exclude_file: t.Optional[Path],
#     sample_metadata_xlsx: Path
# ) -> t.DefaultDict[str, t.List[Path]]:
#     """
#     Processes normal BAM files by validating them, optionally excluding samples present in the exclude file,
#     and splitting the BAM list based on sample sex using the metadata Excel spreadsheet.

#     Args:
#         normal_bams (List[Path]): List of paths to normal BAM files.
#         exclude_file (Optional[Path]): Path to a file specifying samples to exclude.
#         sample_metadata_xlsx (Path): Path to the Excel file containing sample metadata.

#     Returns:
#         dict: A dictionary categorizing BAM files by sample sex.
#     """

#     # Validate the provided BAM files
#     validated_bams = validate_bam_files(normal_bams)
#     logging.debug(f"Validated BAM files: {validated_bams}")

#     # Determine whether to filter BAMs based on the presence of an exclude file
#     if exclude_file:
#         logging.info('Exclude file detected. Filtering normal BAMs...')
#         filtered_bams = remove_unwanted_sample_files(validated_bams, exclude_file)
#         logging.debug(f"Filtered BAM files: {filtered_bams}")
#     else:
#         logging.info('No exclude file detected. Using unfiltered normal BAMs...')
#         filtered_bams = validated_bams

#     # Split the BAM files by sample sex using the provided metadata
#     split_bams = split_file_list_by_sample_sex(filtered_bams, sample_metadata_xlsx)
#     logging.debug(f"Split BAM files by sex: {split_bams}")

#     return split_bams


def main():
    # Get command line arguments
    logging.info("Getting command line arguments...")
    args = parse_arguments()

    # normal_bams = args.normal_bams
    # exclude_file = args.exclude_file
    config_file = args.config
    outdir = args.outdir

    # logging.debug(f"Normal BAMs: {normal_bams}")
    # logging.debug(f"Exclude file: {exclude_file}")
    logging.debug(f"Config file: {config_file}")
    logging.debug(f"Outdir: {outdir}")

    # Extract metadata files from config file
    metadata = extract_metadata_files_from_config_json(config_file)
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

    # Process normal BAMs to get a dictionary of BAMs seperated by sex
    # sex_seperated_normal_bams = process_normal_bams(normal_bams, exclude_file, sample_metadata_xlsx)
    sex_seperated_normal_bams = split_file_list_by_sample_sex(
        normal_bams, sample_metadata_xlsx
    )
    logging.debug(f"Split normal BAM files by sex: {sex_seperated_normal_bams}")

    # Use CNVKit to generate a copy number reference for each sex
    for sex in sex_seperated_normal_bams:
        sex_outdir = Path(outdir / sex)
        sex_outdir.mkdir(parents=True, exist_ok=True)

        # Extract the sex-seperated normal BAMs for the given sex from the dictionary
        sex_normal_bams = sex_seperated_normal_bams[sex]

        # Run CNVKit coverage on the sex separated BAMs
        coverage_files = []
        for bam in sex_normal_bams:
            target_coverage_file = run_cnvkit_coverage(bam, targets_bed, sex_outdir)
            antitarget_coverage_file = run_cnvkit_coverage(
                bam, antitargets_bed, sex_outdir
            )
            coverage_files += [target_coverage_file, antitarget_coverage_file]

        logging.debug(f"Coverage files: {coverage_files}")

        # Run CNVKit reference to create a new copy number reference using the coverage files


if __name__ == "__main__":
    configure_logging()
    main()
