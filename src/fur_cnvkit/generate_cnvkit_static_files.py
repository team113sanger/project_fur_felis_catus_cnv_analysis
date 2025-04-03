import argparse
import json
from pathlib import Path
import typing as t

from fur_cnvkit import constants
from fur_cnvkit.utils.cnvkit_utils import run_cnvkit_access, run_cnvkit_autobin
from fur_cnvkit.utils.file_format_checker import is_fasta, is_bed, validate_bam_files
from fur_cnvkit.utils.fur_utils import (
    remove_unwanted_sample_files,
    categorise_files_by_tumour_normal_status,
)

from fur_cnvkit.utils.logging_utils import setup_logging, get_package_logger

COMMAND_NAME: str = constants.COMMAND_NAME__GENERATE_STATIC_FILES

# Set up the logger
logger = get_package_logger()


def get_argparser(
    subparser: t.Optional[argparse._SubParsersAction] = None,
) -> argparse.ArgumentParser:
    """
    Either returns a new ArgumentParser instance or a subparser for the
    generate_cnvkit_static_files command.

    It is preferrable to use the subparser argument as it unifies the CLI to a
    single entrypoint. To preserve backwards compatibility, the function can
    also be called without the subparser argument.
    """
    if subparser is None:
        parser = argparse.ArgumentParser(
            description=constants.DESCRIPTION__GENERATE_STATIC_FILES
        )
    else:
        parser = subparser.add_parser(
            COMMAND_NAME,
            description=constants.DESCRIPTION__GENERATE_STATIC_FILES,
            help=constants.SHORT_HELP__GENERATE_STATIC_FILES,
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
        "-u",
        metavar="UNPLACED_CONTIGS",
        type=str,
        nargs="*",
        help="Prefixes of unplaced contigs in the reference genome FASTA file. Please provide separated by spaces.",
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

    return parser


def generate_baitset_genes_file(targets_bed: Path, outdir: Path) -> Path:
    logger.info("Generating baitset genes file ...")

    # Extract a set of baitset genes from the targets BED file
    with open(targets_bed, "r") as f:
        baitset_genes = {line.split("\t")[3].strip() for line in f.readlines()}

    # Unpack gene names that have been merged during target BED generation
    unpacked_genes_set = set()

    for gene_string in baitset_genes:
        if "," in gene_string:
            unpacked_genes_set.update(gene_string.split(","))
        else:
            unpacked_genes_set.add(gene_string)

    logger.debug(
        f"Detected {len(unpacked_genes_set)} unique genes : {sorted(unpacked_genes_set)}"
    )

    # Define output file path
    baitset_genes_file_name = targets_bed.name.replace(
        ".target.bed", ".baitset_genes.txt"
    )
    baitset_genes_file_path = outdir / baitset_genes_file_name

    # Write the baitset genes to a file
    with open(baitset_genes_file_path, "w") as f:
        for gene in sorted(unpacked_genes_set):
            f.write(gene + "\n")

    logger.info(f"Baitset genes file saved to {str(baitset_genes_file_path)}")

    return baitset_genes_file_path


def generate_parameter_file(
    bam_files: t.List[Path],
    reference_fasta: Path,
    unplaced_contig_prefixes: t.List[str],
    baitset_bed: Path,
    refflat_file: Path,
    sample_metadata_xlsx: Path,
    access_bed: Path,
    targets_bed: Path,
    antitargets_bed: Path,
    baitset_genes_file: Path,
    parameter_file_name: str,
    outdir: Path,
) -> Path:
    logger.info("Generating parameter file containing CNVKit static files...")

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
        "unplaced_contig_prefixes": unplaced_contig_prefixes,
        "baitset_bed": str(baitset_bed),
        "refflat_file": str(refflat_file),
        "sample_metadata_xlsx": str(sample_metadata_xlsx),
        "access_bed": str(access_bed),
        "targets_bed": str(targets_bed),
        "antitargets_bed": str(antitargets_bed),
        "baitset_genes_file": str(baitset_genes_file),
    }

    # Construct output parameter file path
    output_parameter_file = outdir / f"{parameter_file_name}.parameters.json"

    # Write data to output parameter file path
    logger.info(f"Writing parameter file to {str(output_parameter_file)}")
    try:
        with output_parameter_file.open("w") as json_file:
            json.dump(parameter_data, json_file, indent=4)
            json_file.write("\n")
        logger.info("Successfully wrote parameter file.")
    except Exception as e:
        logger.warning(f"Error writing parameter file: {e}")
        raise

    return output_parameter_file


def main(args: t.Optional[argparse.Namespace] = None) -> None:
    # Get command line arguments
    logger.info("Starting generation of CNVKit static files ...")
    logger.info("Getting and validating command line arguments...")
    if args is None:
        argparser = get_argparser()
        args = argparser.parse_args()
    logger.debug(f"Parsed arguments: {args}")

    validated_bams = validate_bam_files(args.b)
    exclude_file = args.e
    if is_fasta(args.f):
        reference_fasta = args.f
    else:
        raise ValueError(f"{str(args.f)} is not a valid FASTA file.")
    unplaced_contig_prefixes = args.u
    if is_bed(args.t):
        baitset_bed = args.t
    else:
        raise ValueError(f"{str(args.t)} is not a valid FASTA file.")
    refflat_file = args.r
    sample_metadata_xlsx = args.m
    parameter_file_name = args.p
    outdir = args.o

    logger.debug(f"BAM files: {validated_bams}")
    logger.debug(f"Reference FASTA: {reference_fasta}")
    logger.debug(f"Baitset BED: {baitset_bed}")
    logger.debug(f"RefFlat file: {refflat_file}")
    logger.debug(f"Sample metadata Excel spreadsheet: {sample_metadata_xlsx}")
    logger.debug(f"Outdir: {outdir}")

    logger.info(
        "Command line arguments will now be used to generate CNVKit statis files ..."
    )

    # Remove any unwanted samples from final BAM list
    if exclude_file:
        logger.info("Exclude file detected. Filtering normal BAMs...")
        filtered_bams = remove_unwanted_sample_files(validated_bams, exclude_file)
        logger.debug(f"Filtered BAM files: {filtered_bams}")
    else:
        logger.info("No exclude file detected. Using unfiltered normal BAMs...")
        filtered_bams = validated_bams

    # Run cnvkit.py access
    access_bed = run_cnvkit_access(reference_fasta, outdir)

    # Run cnvkit.py autobin
    target_bed_dict = run_cnvkit_autobin(
        filtered_bams, baitset_bed, access_bed, refflat_file, outdir
    )

    # Generate a baitset genes file using the target BED
    baitset_genes_file = generate_baitset_genes_file(target_bed_dict["target"], outdir)

    # Generate parameter file including newly generated reference files
    generate_parameter_file(
        filtered_bams,
        reference_fasta,
        unplaced_contig_prefixes,
        baitset_bed,
        refflat_file,
        sample_metadata_xlsx,
        access_bed,
        target_bed_dict["target"],
        target_bed_dict["antitarget"],
        baitset_genes_file,
        parameter_file_name,
        outdir,
    )

    logger.info("CNVKit static file generation successfully completed.")


if __name__ == "__main__":
    setup_logging()
    main()
