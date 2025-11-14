import typing as t
import argparse
import os
from pathlib import Path
import multiprocessing

from fur_cnvkit import constants
from fur_cnvkit.normal_vs_normal import perform_normal_vs_normal_comparisons
from fur_cnvkit.utils.cnvkit_utils import (
    run_cnvkit_coverage,
    run_cnvkit_reference,
    run_cnvkit_sex,
)
from fur_cnvkit.utils.fur_utils import (
    extract_metadata_files_from_parameter_json,
    split_file_list_by_sample_sex,
    get_sample_ids_for_file_list,
    set_metadata_columns,
)
from fur_cnvkit.utils.logging_utils import setup_logging, get_package_logger

# Constants
COMMAND_NAME: str = constants.COMMAND_NAME__GENERATE_CN_REFERENCE

# Set up logging
logger = get_package_logger()


def get_argparser(  # noqa: C901
    subparser: t.Optional[argparse._SubParsersAction] = None,
) -> argparse.ArgumentParser:
    """
    Either returns a new ArgumentParser instance or a subparser for the
    generate_copy_number_reference command.

    It is preferrable to use the subparser argument as it unifies the CLI to a
    single entrypoint. To preserve backwards compatibility, the function can
    also be called without the subparser argument.
    """
    if subparser is None:
        parser = argparse.ArgumentParser(
            description=constants.DESCRIPTION__GENERATE_CN_REFERENCE
        )
    else:
        parser = subparser.add_parser(
            COMMAND_NAME,
            description=constants.DESCRIPTION__GENERATE_CN_REFERENCE,
            help=constants.SHORT_HELP__GENERATE_CN_REFERENCE,
        )

    def _parse_max_cpus(value: t.Union[int, str]) -> int:
        """
        Parse and validate the max_cpus argument.

        Args:
            value: The value provided by the user, or None if --max_cpus was specified without a value

        Returns:
            int: The number of CPUs to use

        Raises:
            argparse.ArgumentTypeError: If the provided value is invalid
        """
        max_available = multiprocessing.cpu_count()

        try:
            cpu_count = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid value for --max_cpus: {value}. Must be an integer."
            )

        if cpu_count <= 0:
            raise argparse.ArgumentTypeError(
                f"Invalid value for --max_cpus: {value}. Must be a positive integer."
            )
        if cpu_count > max_available:
            raise argparse.ArgumentTypeError(
                f"Invalid value for --max_cpus: {value}. Must be less than {max_available}."
            )

        return cpu_count

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
    parser.add_argument(
        "--verbose",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        dest="verbose",
        help="log level",
    )
    parser.add_argument(
        "--max_cpus",
        type=_parse_max_cpus,
        nargs="?",  # Makes the argument optional after the flag
        const=multiprocessing.cpu_count(),  # This value is used when --max_cpus is specified without a value
        default=1,  # This value is used when --max_cpus is not specified at all
        dest="max_cpus",
        help="Maximum number of CPUs to use. If specified without a value, all available CPUs will be used. Default is 1.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Maximum number of normal-vs-normal worker processes to run concurrently. Defaults to min(requested CPUs, available CPUs).",
    )
    return parser


def reclassify_unknown_samples(
    sex_separated_normal_bams: dict,
    targets_bed: Path,
    antitargets_bed: Path,
    outdir: Path,
    max_cpus: t.Optional[int] = None,
) -> dict:
    """
    For any samples initially marked as 'unknown' (from metadata), generate
    coverage files and use run_cnvkit_sex to determine their sex. Each sample's
    assignment is recorded in a file in the unknown directory, and the sample
    is added to the appropriate 'male' or 'female' group.

    Returns:
        dict: The updated dictionary grouping samples by sex.
    """
    if "unknown" not in sex_separated_normal_bams:
        return sex_separated_normal_bams

    unknown_samples = sex_separated_normal_bams.pop("unknown")
    unknown_dir = outdir / "unknown"
    unknown_cov_dir = unknown_dir / "coverage_files"
    unknown_cov_dir.mkdir(parents=True, exist_ok=True)

    # Dictionary to store assigned sexes for unknown samples.
    unknown_assignments = {}

    for sample_bam in unknown_samples:
        logger.info(f"Reclassifying sample {sample_bam} with unknown sex ...")
        target_cov = run_cnvkit_coverage(
            sample_bam, targets_bed, unknown_cov_dir, cpus=max_cpus
        )
        antitarget_cov = run_cnvkit_coverage(
            sample_bam, antitargets_bed, unknown_cov_dir, cpus=max_cpus
        )
        determined_sex = run_cnvkit_sex([target_cov, antitarget_cov])
        logger.info(f"Determined sex for sample {sample_bam}: {determined_sex}")
        unknown_assignments[sample_bam.name] = determined_sex
        if determined_sex in ("male", "female"):
            sex_separated_normal_bams.setdefault(determined_sex, []).append(sample_bam)
        else:
            logger.warning(
                f"Sample {sample_bam} remains unknown after reclassification."
            )

    # Write the assignments to a file in the unknown directory.
    assignments_file = unknown_dir / "unknown_assigned_sexes.txt"
    with assignments_file.open("w") as f:
        for sample_name, assigned_sex in unknown_assignments.items():
            f.write(f"{sample_name}\t{assigned_sex}\n")
    logger.info(f"Unknown sample assignments written to {assignments_file}")

    return sex_separated_normal_bams


def generate_reference_for_sex(
    sex: str,
    sex_normal_bams: list,
    outdir: Path,
    reference_fasta: Path,
    targets_bed: Path,
    antitargets_bed: Path,
    sample_metadata_xlsx: Path,
    max_cpus: t.Optional[int] = None,
    threads: t.Optional[int] = None,
):
    """
    For a given sex, generate coverage files, perform normal vs. normal comparisons,
    and create a copy number reference using the filtered coverage files.
    """
    logger.info(f"Running CNVKit reference generation for {sex} samples ...")
    sex_outdir = outdir / sex
    sex_outdir.mkdir(parents=True, exist_ok=True)
    coverage_file_dir = sex_outdir / "coverage_files"
    coverage_file_dir.mkdir(parents=True, exist_ok=True)

    # Generate coverage files for each sample.
    sex_normal_coverage_files = []
    for sample_bam in sex_normal_bams:
        logger.info(f"Generating coverage files for sample {sample_bam} ({sex})...")
        target_cov = run_cnvkit_coverage(
            sample_bam, targets_bed, coverage_file_dir, cpus=max_cpus
        )
        antitarget_cov = run_cnvkit_coverage(
            sample_bam, antitargets_bed, coverage_file_dir, cpus=max_cpus
        )
        sex_normal_coverage_files.extend([target_cov, antitarget_cov])
    logger.debug(f"Coverage files for {sex} samples: {sex_normal_coverage_files}")

    # Perform normal vs. normal comparisons to filter the coverage files.
    normal_vs_normal_dir = sex_outdir / "normal_vs_normal"
    normal_vs_normal_dir.mkdir(parents=True, exist_ok=True)
    filtered_coverage_files = perform_normal_vs_normal_comparisons(
        sex_normal_coverage_files,
        reference_fasta,
        sample_metadata_xlsx,
        normal_vs_normal_dir,
        max_workers=threads,
    )
    logger.debug(f"Filtered coverage files for {sex}: {filtered_coverage_files}")

    # Record the samples used for the reference.
    used_samples_file_path = sex_outdir / "samples_used_in_reference.txt"
    samples_used_in_reference_list = get_sample_ids_for_file_list(
        filtered_coverage_files
    )
    with open(used_samples_file_path, "w") as f:
        for sample_id in samples_used_in_reference_list:
            f.write(f"{sample_id}\n")

    # Generate the copy number reference.
    logger.info(f"Generating copy number reference for {sex} samples ...")
    sex_reference_file = run_cnvkit_reference(
        coverage_files=filtered_coverage_files,
        reference_fasta=reference_fasta,
        output_prefix=sex,
        outdir=sex_outdir,
        sex=sex,
    )
    logger.debug(f"Reference for {sex} samples generated at: {sex_reference_file}")


def main(args: t.Optional[argparse.Namespace] = None):
    if args is None:
        argparser = get_argparser()
        args = argparser.parse_args()
    logger.debug(f"Parsed arguments: {args}")

    parameter_file = args.parameter_file
    outdir = args.outdir
    max_cpus = args.max_cpus
    threads = args.threads
    if threads is not None:
        threads = max(1, min(threads, os.cpu_count() or 1))

    logger.debug(f"Parameter file: {parameter_file}")
    logger.debug(f"Outdir: {outdir}")

    # Extract metadata from the parameter file.
    metadata = extract_metadata_files_from_parameter_json(parameter_file)
    set_metadata_columns(metadata.get("metadata_columns"))
    normal_bams = metadata["normal_bams"]
    reference_fasta = metadata["reference_fasta"]
    baitset_bed = metadata["baitset_bed"]
    sample_metadata_xlsx = metadata["sample_metadata_xlsx"]
    targets_bed = metadata["targets_bed"]
    antitargets_bed = metadata["antitargets_bed"]

    logger.debug(f"Normal BAMs: {normal_bams}")
    logger.debug(f"Reference FASTA: {reference_fasta}")
    logger.debug(f"Baitset BED: {baitset_bed}")
    logger.debug(f"Sample metadata Excel: {sample_metadata_xlsx}")
    logger.debug(f"Targets BED: {targets_bed}")
    logger.debug(f"Antitargets BED: {antitargets_bed}")

    # Group normal BAMs by their metadata sex.
    logger.info("Separating normal BAMs based on their metadata sex ...")
    sex_separated_normal_bams = split_file_list_by_sample_sex(
        normal_bams, sample_metadata_xlsx
    )
    logger.debug(f"Initial grouping by sex: {sex_separated_normal_bams}")

    # Reclassify unknown samples using coverage files.
    sex_separated_normal_bams = reclassify_unknown_samples(
        sex_separated_normal_bams=sex_separated_normal_bams,
        targets_bed=targets_bed,
        antitargets_bed=antitargets_bed,
        outdir=outdir,
        max_cpus=max_cpus,
    )
    logger.debug(f"Final grouping of samples by sex: {sex_separated_normal_bams}")

    # For each sex, generate one copy number reference.
    for sex in ("male", "female"):
        if sex not in sex_separated_normal_bams:
            logger.info(f"No {sex} samples available, skipping reference generation.")
            continue
        generate_reference_for_sex(
            sex=sex,
            sex_normal_bams=sex_separated_normal_bams[sex],
            outdir=outdir,
            reference_fasta=reference_fasta,
            targets_bed=targets_bed,
            antitargets_bed=antitargets_bed,
            sample_metadata_xlsx=sample_metadata_xlsx,
            max_cpus=max_cpus,
            threads=threads,
        )


if __name__ == "__main__":
    setup_logging()
    main()
