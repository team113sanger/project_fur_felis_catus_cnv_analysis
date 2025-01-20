import argparse
import json
import logging
from pathlib import Path
import typing as t

from fur_cnvkit.utils.cnvkit_utils import run_cnvkit_access, run_cnvkit_autobin
from fur_cnvkit.utils.file_format_checker import is_fasta, is_bed, validate_bam_files
from fur_cnvkit.utils.fur_utils import (
    remove_unwanted_sample_files,
    categorise_files_by_tumour_normal_status,
)


def configure_logging():
    """
    Define logging configuration
    """
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate various CNVKit files needed to run downstream analyses."
    )
    parser.add_argument(
        "-b",
        metavar="BAMS",
        type=Path,
        nargs="+",
        required=True,
        help="Path to tumour and normal BAM files.",
    )
    parser.add_argument(
        "-e",
        metavar="EXCLUDE_FILE",
        type=Path,
        required=False,
        help="Path to exclude file containing samples to remove from final analyses",
    )
    parser.add_argument(
        "-f",
        metavar="FASTA",
        type=Path,
        required=True,
        help="Path to reference genome FASTA file.",
    )
    parser.add_argument(
        "-t",
        metavar="BAITSET",
        type=Path,
        required=True,
        help="Path to baitset BED file.",
    )
    parser.add_argument(
        "-r",
        metavar="REFFLAT",
        type=Path,
        required=True,
        help="Path to refFlat annotation file.",
    )
    parser.add_argument(
        "-m",
        metavar="METADATA",
        type=Path,
        required=True,
        help="Path to sample metadata Excel spreadsheet.",
    )
    parser.add_argument(
        "-p",
        metavar="PARAMETER_FILE_PREFIX",
        type=str,
        required=False,
        default="parameters.json",
        help="Prefix of output parameter file to store CNVKit static files ({prefix}.parameters.json). Can be used in downstream scripts to avoid re-specifying paths. Will be written to outdir (-o). Default: parameters.json",
    )
    parser.add_argument(
        "-o",
        metavar="OUTDIR",
        type=Path,
        required=True,
        help="Path to directory where output files will be stored.",
    )

    return parser.parse_args()


def generate_parameter_file(
    bam_files: t.List[Path],
    reference_fasta: Path,
    baitset_bed: Path,
    refflat_file: Path,
    sample_metadata_xlsx: Path,
    targets_bed: Path,
    antitargets_bed: Path,
    parameter_file_name: str,
    outdir: Path,
) -> Path:
    logging.info("Generating parameter file containing CNVKit static files...")

    # Categorise the BAM files based on their tumour/normal status
    tn_status_bam_dict = categorise_files_by_tumour_normal_status(
        bam_files, sample_metadata_xlsx
    )
    tumour_bams = [str(bam) for bam in tn_status_bam_dict["T"]]
    normal_bams = [str(bam) for bam in tn_status_bam_dict["N"]]

    # Initialise a dictionary containing the parameter file data
    parameter_data = {
        "all_bams": [str(bam) for bam in bam_files],
        "tumour_bams": tumour_bams,
        "normal_bams": normal_bams,
        "reference_fasta": str(reference_fasta),
        "baitset_bed": str(baitset_bed),
        "refflat_file": str(refflat_file),
        "sample_metadata_xlsx": str(sample_metadata_xlsx),
        "targets_bed": str(targets_bed),
        "antitargets_bed": str(antitargets_bed),
    }

    # Construct output parameter file path
    output_parameter_file = outdir / f"{parameter_file_name}.parameters.json"

    # Write data to output parameter file path
    logging.info(f"Writing parameter file to {str(output_parameter_file)}")
    try:
        with output_parameter_file.open("w") as json_file:
            json.dump(parameter_data, json_file, indent=4)
            json_file.write("\n")
        logging.info("Successfully wrote parameter file.")
    except Exception as e:
        logging.warning(f"Error writing parameter file: {e}")
        raise

    return output_parameter_file


def main():
    # Get command line arguments
    logging.info("Starting generation of CNVKit static files ...")
    logging.info("Getting and validating command line arguments...")
    args = parse_arguments()

    validated_bams = validate_bam_files(args.b)
    exclude_file = args.e
    if is_fasta(args.f):
        reference_fasta = args.f
    else:
        raise ValueError(f"{str(args.f)} is not a valid FASTA file.")
    if is_bed(args.t):
        baitset_bed = args.t
    else:
        raise ValueError(f"{str(args.t)} is not a valid FASTA file.")
    refflat_file = args.r
    sample_metadata_xlsx = args.m
    parameter_file_name = args.p
    outdir = args.o

    logging.debug(f"BAM files: {validated_bams}")
    logging.debug(f"Reference FASTA: {reference_fasta}")
    logging.debug(f"Baitset BED: {baitset_bed}")
    logging.debug(f"RefFlat file: {refflat_file}")
    logging.debug(f"Sample metadata Excel spreadsheet: {sample_metadata_xlsx}")
    logging.debug(f"Outdir: {outdir}")

    logging.info(
        "Command line arguments will now be used to generate CNVKit statis files ..."
    )

    # Remove any unwanted samples from final BAM list
    if exclude_file:
        logging.info("Exclude file detected. Filtering normal BAMs...")
        filtered_bams = remove_unwanted_sample_files(validated_bams, exclude_file)
        logging.debug(f"Filtered BAM files: {filtered_bams}")
    else:
        logging.info("No exclude file detected. Using unfiltered normal BAMs...")
        filtered_bams = validated_bams

    # Run cnvkit.py access
    access_bed = run_cnvkit_access(reference_fasta, outdir)

    # Run cnvkit.py autobin
    target_bed_dict = run_cnvkit_autobin(
        filtered_bams, baitset_bed, access_bed, refflat_file, outdir
    )

    # Generate parameter file including newly generated reference files
    generate_parameter_file(
        filtered_bams,
        reference_fasta,
        baitset_bed,
        refflat_file,
        sample_metadata_xlsx,
        target_bed_dict["target"],
        target_bed_dict["antitarget"],
        parameter_file_name,
        outdir,
    )

    logging.info("CNVKit static file generation successfully completed.")


if __name__ == "__main__":
    configure_logging()
    main()
