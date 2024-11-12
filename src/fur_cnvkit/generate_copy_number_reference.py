import argparse
import logging
from pathlib import Path

from utils.file_format_checker import validate_bam_files


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
        "-n",
        "--normal_bams",
        type=Path,
        nargs="+",
        required=True,
        help="Path to normal BAM(s) to generate copy number referemce with.",
    )
    parser.add_argument(
        "-e",
        "--exclude_file",
        type=Path,
        nargs=1,
        required=False,
        help="Exclude file containing samples to exclude from copy number reference (One sample on each line).",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        nargs=1,
        required=True,
        help="Path to the configuration file.",
    )

    return parser.parse_args()


def main():
    # Get command line arguments
    logging.info("Getting command line arguments...")
    args = parse_arguments()

    normal_bams = validate_bam_files(args.normal_bams)
    exclude_file = args.exclude_file[0]
    config_file = args.config_file[0]

    logging.debug(f"Normal BAMs: {normal_bams}")
    logging.debug(f"Exclude file: {exclude_file}")
    logging.debug(f"Config file: {config_file}")

    # Step 1 - Create a


if __name__ == "__main__":
    configure_logging()
    main()
