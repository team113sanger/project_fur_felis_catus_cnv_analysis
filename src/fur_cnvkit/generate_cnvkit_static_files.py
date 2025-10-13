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
    set_metadata_columns,
    DEFAULT_METADATA_COLUMNS,
)

from fur_cnvkit.utils.logging_utils import setup_logging, get_package_logger

COMMAND_NAME: str = constants.COMMAND_NAME__GENERATE_STATIC_FILES
DEFAULT_PARAMETER_FILE_NAME: str = "parameters.json"

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
        help="Path to sample metadata file (Excel, TSV, or CSV).",
    )
    parser.add_argument(
        "--metadata-sample-id-column",
        type=str,
        default=DEFAULT_METADATA_COLUMNS.sample_id,
        help=(
            "Column name containing sample identifiers in the metadata file. "
            f"Default: '{DEFAULT_METADATA_COLUMNS.sample_id}'."
        ),
    )
    parser.add_argument(
        "--metadata-tumour-normal-column",
        type=str,
        default=DEFAULT_METADATA_COLUMNS.tumour_normal,
        help=(
            "Column name describing tumour/normal status. "
            f"Default: '{DEFAULT_METADATA_COLUMNS.tumour_normal}'."
        ),
    )
    parser.add_argument(
        "--metadata-sex-column",
        type=str,
        default=DEFAULT_METADATA_COLUMNS.sex,
        help=(
            "Column name describing sample sex (values should be M/F/U). "
            f"Default: '{DEFAULT_METADATA_COLUMNS.sex}'."
        ),
    )
    parser.add_argument(
        "--metadata-study-column",
        type=str,
        default=DEFAULT_METADATA_COLUMNS.study,
        help=(
            "Optional column name describing study/group identifiers. "
            "If omitted, Excel sheet names or a single 'default' cohort will be used."
        ),
    )
    parser.add_argument(
        "-p",
        metavar="PARAMETER_FILE_PREFIX",
        type=str,
        required=False,
        default=DEFAULT_PARAMETER_FILE_NAME,
        help=(
            f"Prefix of output parameter file to store CNVKit static files ({{prefix}}.parameters.json). "
            "Can be used in downstream scripts to avoid re-specifying paths. "
            f"Will be written to outdir (-o). Default: {DEFAULT_PARAMETER_FILE_NAME}"
        ),
    )
    parser.add_argument(
        "-x",
        metavar="ACCESS_EXCLUDE_BED",
        type=Path,
        required=False,
        help="BED file of regions to exclude when computing accessible regions "
        "(passed to `cnvkit.py access -x`).",
    )
    parser.add_argument(
        "-o",
        metavar="OUTDIR",
        type=Path,
        required=True,
        help="Path to directory where output files will be stored.",
    )
    parser.add_argument(
        "--verbose",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        dest="verbose",
        help="log level",
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
    metadata_columns: t.Optional[t.Dict[str, t.Optional[str]]] = None,
) -> Path:
    logger.info("Generating parameter file containing CNVKit static files...")

    if metadata_columns is not None:
        set_metadata_columns(metadata_columns)

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
    if metadata_columns:
        metadata_columns_filtered = {
            key: value for key, value in metadata_columns.items() if value is not None
        }
        if metadata_columns_filtered:
            parameter_data["metadata_columns"] = metadata_columns_filtered

    # Construct output parameter file path (default: parameters.json from argparse)
    endswith_default = parameter_file_name.endswith(DEFAULT_PARAMETER_FILE_NAME)
    is_default = parameter_file_name == DEFAULT_PARAMETER_FILE_NAME
    if endswith_default and not is_default:
        # E.g "my_parameter_file.parameters.json" -> "my_parameter_file.parameters.json"
        prefix = parameter_file_name.replace(DEFAULT_PARAMETER_FILE_NAME, "").rstrip(
            "."
        )
        parameter_file_name = f"{prefix}.{DEFAULT_PARAMETER_FILE_NAME}"
    elif not endswith_default:
        # E.g "my_parameter_file" -> "my_parameter_file.parameters.json"
        prefix = parameter_file_name.rstrip(".")
        parameter_file_name = f"{prefix}.{DEFAULT_PARAMETER_FILE_NAME}"

    output_parameter_file = outdir / parameter_file_name

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


def validate_fasta(fasta: Path) -> Path:
    if not is_fasta(fasta):
        raise ValueError(f"{str(fasta)} is not a valid FASTA file. Check input data.")
    return fasta


def validate_bed(bed: Path) -> Path:
    if not is_bed(bed):
        raise ValueError(f"{str(bed)} is not a valid BED file. Check input data.")
    return bed


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
    reference_fasta = validate_fasta(args.f)
    unplaced_contig_prefixes = args.u
    baitset_bed = validate_bed(args.t)
    refflat_file = args.r
    sample_metadata_xlsx = args.m
    parameter_file_name = args.p
    access_exclude_bed = validate_bed(args.x) if args.x else None
    outdir = args.o
    metadata_columns_input = {
        "sample_id": args.metadata_sample_id_column,
        "tumour_normal": args.metadata_tumour_normal_column,
        "sex": args.metadata_sex_column,
        "study": args.metadata_study_column,
    }
    # Configure metadata parsing for downstream utilities.
    set_metadata_columns(metadata_columns_input)
    default_metadata_columns = {
        "sample_id": DEFAULT_METADATA_COLUMNS.sample_id,
        "tumour_normal": DEFAULT_METADATA_COLUMNS.tumour_normal,
        "sex": DEFAULT_METADATA_COLUMNS.sex,
        "study": DEFAULT_METADATA_COLUMNS.study,
    }
    metadata_columns_clean = {
        key: value for key, value in metadata_columns_input.items() if value is not None
    }
    default_columns_clean = {
        key: value
        for key, value in default_metadata_columns.items()
        if value is not None
    }
    metadata_columns_for_parameter = (
        metadata_columns_clean
        if metadata_columns_clean != default_columns_clean
        else None
    )

    logger.debug(f"BAM files: {validated_bams}")
    logger.debug(f"Reference FASTA: {reference_fasta}")
    logger.debug(f"Baitset BED: {baitset_bed}")
    logger.debug(f"RefFlat file: {refflat_file}")
    logger.debug(f"Sample metadata file: {sample_metadata_xlsx}")
    logger.debug(f"Metadata column configuration: {metadata_columns_input}")
    logger.debug(f"Parameter file name: {parameter_file_name}")
    if access_exclude_bed is not None:
        logger.debug(f"Access exclude BED: {access_exclude_bed}")
    else:
        logger.debug("No access exclude BED provided.")
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
    access_bed = run_cnvkit_access(
        reference_fasta, outdir, exclude_bed=access_exclude_bed
    )

    # Run cnvkit.py autobin
    target_bed_dict = run_cnvkit_autobin(
        filtered_bams, baitset_bed, access_bed, refflat_file, outdir
    )

    # Generate a baitset genes file using the target BED
    baitset_genes_file = generate_baitset_genes_file(target_bed_dict["target"], outdir)

    # Generate parameter file including newly generated reference files
    generate_parameter_file(
        bam_files=filtered_bams,
        reference_fasta=reference_fasta,
        unplaced_contig_prefixes=unplaced_contig_prefixes,
        baitset_bed=baitset_bed,
        refflat_file=refflat_file,
        sample_metadata_xlsx=sample_metadata_xlsx,
        access_bed=access_bed,
        targets_bed=target_bed_dict["target"],
        antitargets_bed=target_bed_dict["antitarget"],
        baitset_genes_file=baitset_genes_file,
        metadata_columns=metadata_columns_for_parameter,
        parameter_file_name=parameter_file_name,
        outdir=outdir,
    )

    logger.info("CNVKit static file generation successfully completed.")


if __name__ == "__main__":
    setup_logging()
    main()
