import argparse
from pathlib import Path

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
)
from fur_cnvkit.utils.logging_utils import setup_logging, get_package_logger

# Set up logging
logger = get_package_logger()


def parse_arguments():
    """Define and parse command line arguments."""
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


def reclassify_unknown_samples(
    sex_separated_normal_bams: dict,
    targets_bed: Path,
    antitargets_bed: Path,
    outdir: Path,
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
        target_cov = run_cnvkit_coverage(sample_bam, targets_bed, unknown_cov_dir)
        antitarget_cov = run_cnvkit_coverage(
            sample_bam, antitargets_bed, unknown_cov_dir
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
        target_cov = run_cnvkit_coverage(sample_bam, targets_bed, coverage_file_dir)
        antitarget_cov = run_cnvkit_coverage(
            sample_bam, antitargets_bed, coverage_file_dir
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


def main():
    logger.info("Getting command line arguments...")
    args = parse_arguments()
    parameter_file = args.parameter_file
    outdir = args.outdir

    logger.debug(f"Parameter file: {parameter_file}")
    logger.debug(f"Outdir: {outdir}")

    # Extract metadata from the parameter file.
    metadata = extract_metadata_files_from_parameter_json(parameter_file)
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
        sex_separated_normal_bams, targets_bed, antitargets_bed, outdir
    )
    logger.debug(f"Final grouping of samples by sex: {sex_separated_normal_bams}")

    # For each sex, generate one copy number reference.
    for sex in ("male", "female"):
        if sex not in sex_separated_normal_bams:
            logger.info(f"No {sex} samples available, skipping reference generation.")
            continue
        generate_reference_for_sex(
            sex,
            sex_separated_normal_bams[sex],
            outdir,
            reference_fasta,
            targets_bed,
            antitargets_bed,
            sample_metadata_xlsx,
        )


if __name__ == "__main__":
    setup_logging()
    main()
