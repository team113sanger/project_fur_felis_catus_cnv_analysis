import logging
from pathlib import Path
import re
import shlex
import subprocess
import typing as t
import os

import pandas as pd

from fur_cnvkit.utils.file_format_checker import is_bam, is_bed, is_fasta

logger = logging.getLogger(__name__)


# Helper functions for command execution
def run_command(command: str, env: t.Optional[dict[str, str]] = None) -> None:
    """Run a given command using subprocess, including logging statements and exception handling.

    Args:
        command (str): Command to be executed.
    """
    logger.debug(f"Executing command: {command}")
    try:
        result = execute_command(command, env)
        log_success(command, result)
    except subprocess.CalledProcessError as e:
        log_error(command, e)


def execute_command(
    command: str, env: t.Optional[dict[str, str]] = None
) -> subprocess.CompletedProcess:
    """Execute a command using subprocess and return the result.

    Args:
        command (str): Command to execute.

    Returns:
        subprocess.CompletedProcess: The result of the executed command.

    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    if env is None:
        env = os.environ.copy()
    command_list = shlex.split(command)
    return subprocess.run(
        command_list,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
        env=env,
    )


def log_success(command: str, result: subprocess.CompletedProcess) -> None:
    """Log the success of a command execution.

    Args:
        command (str): The executed command.
        result (subprocess.CompletedProcess): The result of the command execution.
    """
    logger.info(f"Successfully completed command: {command}")
    if result.stdout:
        logger.debug(f"Command Standard Output: {result.stdout.strip()}")
    if result.stderr:
        logger.debug(f"Command Standard Error: {result.stderr.strip()}")


def log_error(command: str, error: subprocess.CalledProcessError) -> None:
    """Log details of a failed command execution.

    Args:
        command (str): The command that failed.
        error (subprocess.CalledProcessError): The exception raised.
    """
    logger.error(f"Command '{command}' failed with exit code {error.returncode}")
    if error.stdout:
        logger.error(f"Error Output: {error.stdout.strip()}")
    if error.stderr:
        logger.error(f"Error Details: {error.stderr.strip()}")


def convert_file_list_to_string(files: t.List[Path]) -> str:
    if len(files) == 0:
        logger.warning("No files detected. Please check input data.")
        return ""
    elif len(files) == 1:
        file_path_string = str(files[0].absolute())
    else:
        file_path_string = " ".join([str(file.absolute()) for file in files])

    return file_path_string


# -----------------------------------------------------------------------------
# Helper function for skipping file generation if output already exists
# -----------------------------------------------------------------------------
def skip_file_generation(
    file: Path,
    validator: t.Optional[t.Callable[[Path], bool]] = None,
    skip_message: t.Optional[str] = None,
) -> bool:
    """
    Check if a file already exists and optionally validate its contents.
    If the file exists and (if provided) passes the validator, log a message and return True
    (indicating that file generation can be skipped). Otherwise, return False.

    Args:
        file (Path): The file to check.
        validator (Optional[Callable[[Path], bool]]): A function that takes the file path and returns True
            if the file's contents are as expected. If None, only file existence is checked.
        skip_message (Optional[str]): Custom log message if the file exists and is valid.

    Returns:
        bool: True if file exists (and is valid when validator is provided), otherwise False.
    """
    if file.exists():
        if validator and not validator(file):
            logger.warning(f"File {file} exists but did not pass validation.")
            return False
        msg = (
            skip_message
            or f"File {file} already exists and is valid. Skipping generation."
        )
        logger.info(msg)
        return True
    return False


# -----------------------------------------------------------------------------
# Validator functions for output files
# -----------------------------------------------------------------------------
def file_exists(file: Path) -> bool:
    """
    Check if the file exists.

    Args:
        file (Path): The file to check.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    if not file.exists():
        logger.warning(
            f"The file '{file}' does not exist. Please check the input path."
        )
        return False
    return True


def file_exists_and_nonempty(path: Path, min_size: int = 1) -> bool:
    """
    Returns True if 'path' exists and has size >= min_size bytes.
    Otherwise False.
    """
    if not path.exists():
        return False
    if path.stat().st_size < min_size:
        return False
    return True


def is_regular_file(file: Path) -> bool:
    """
    Check if the path is a regular file.

    Args:
        file (Path): The file to check.

    Returns:
        bool: True if it's a regular file, False otherwise.
    """
    if not file.is_file():
        logger.warning(
            f"The path '{file}' is not a valid file. Please check the input."
        )
        return False
    return True


def has_correct_suffix(file: Path, expected_suffix: str) -> bool:
    """
    Check if the file has the expected suffix.

    Args:
        file (Path): The file to check.
        expected_suffix (str): The expected file suffix.

    Returns:
        bool: True if the file has the expected suffix, False otherwise.
    """
    if not file.name.endswith(expected_suffix):
        logger.warning(
            f"The file '{file}' has an unexpected suffix. Expected a '{expected_suffix}' file."
        )
        return False
    return True


def has_expected_header(file: Path, expected_header: t.Union[str, t.List[str]]) -> bool:
    """
    Check if the file has the expected header or one of the valid headers.

    Args:
        file (Path): The file to check.
        expected_header (Union[str, List[str]]): The expected header line as a string,
            or a list of valid header lines.

    Returns:
        bool: True if the file has one of the expected headers, False otherwise.
    """
    try:
        with file.open() as f:
            actual_header = f.readline().strip()

        # Normalize expected_header to a list, in case it's a single string.
        valid_headers = (
            expected_header if isinstance(expected_header, list) else [expected_header]
        )

        if actual_header in valid_headers:
            logger.info(f"Successfully validated the structure of '{file}'.")
            return True
        else:
            logger.warning(
                f"The file '{file}' has unexpected columns. Expected one of: {valid_headers} Got: '{actual_header}'"
            )
            return False
    except Exception as e:
        logger.error(f"An error occurred while reading the file '{file}': {e}")
        return False


def is_valid_coverage_file(file: Path) -> bool:
    """
    Check if an input coverage file has the expected structure.
    Expected suffix: "coverage.cnn"
    Expected header: "chromosome\tstart\tend\tgene\tdepth\tlog2"
    """
    return (
        file_exists(file)
        and is_regular_file(file)
        and has_correct_suffix(file, "coverage.cnn")
        and has_expected_header(file, "chromosome\tstart\tend\tgene\tdepth\tlog2")
    )


def is_valid_targetcoverage_file(file: Path) -> bool:
    """Validate that the target coverage file exists, is non-empty, and has the expected suffix."""
    return file_exists_and_nonempty(file) and has_correct_suffix(
        file, ".targetcoverage.cnn"
    )


def is_valid_antitargetcoverage_file(file: Path) -> bool:
    """Validate that the antitarget coverage file exists, is non-empty, and has the expected suffix."""
    return file_exists_and_nonempty(file) and has_correct_suffix(
        file, ".antitargetcoverage.cnn"
    )


def validate_sample_coverage_files(
    sample_id: str, sample_coverage_files: t.List[Path]
) -> t.Tuple[Path, Path]:
    """
    Validates that a sample has exactly one target coverage file and one antitarget coverage file.

    Args:
        sample_id (str): The unique identifier of the sample.
        sample_coverage_files (List[Path]): A list of file paths to check.

    Returns:
        Tuple[Path, Path]: A tuple containing the target coverage file and the antitarget coverage file.

    Raises:
        ValueError: If the expected files are missing or their counts are incorrect.

    Logs:
        Logs warnings for unexpected counts and information about successfully identified files.
    """

    def filter_files(file_list: t.List[Path], keyword: str) -> t.List[Path]:
        """Filter files containing the specified keyword in their names."""
        return [file for file in file_list if keyword in file.name]

    target_files = filter_files(sample_coverage_files, ".targetcoverage.cnn")
    antitarget_files = filter_files(sample_coverage_files, ".antitargetcoverage")

    if not target_files or not antitarget_files:
        raise ValueError(
            f"Insufficient coverage files for sample '{sample_id}'. "
            f"Expected one target and one antitarget file, but found: "
            f"Target files: {len(target_files)}, Antitarget files: {len(antitarget_files)}."
        )

    if len(target_files) > 1 or len(antitarget_files) > 1:
        logger.warning(
            f"Unexpected number of coverage files for sample '{sample_id}'. "
            f"Target files: {target_files}, Antitarget files: {antitarget_files}."
        )

    logger.info(
        f"Successfully identified coverage files for sample '{sample_id}'. "
        f"Target file: {target_files[0]}, Antitarget file: {antitarget_files[0]}."
    )

    return target_files[0], antitarget_files[0]


def is_valid_reference_file(file: Path) -> bool:
    """Validate that the reference file exists, is non-empty, and has a '.reference.cnn' suffix."""
    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and has_correct_suffix(file, ".reference.cnn")
        and has_expected_header(
            file, "chromosome\tstart\tend\tgene\tlog2\tdepth\tgc\trmask\tspread"
        )
    )


def is_valid_cnr_file(file: Path) -> bool:
    """Validate that the CNV ratio file exists, is non-empty, has a '.cnr' suffix, and has the expected header."""
    valid_headers = [
        "chromosome\tstart\tend\tgene\tdepth\tlog2\tweight",
        "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight",
    ]

    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and has_correct_suffix(file, ".cnr")
        and has_expected_header(file, valid_headers)
    )


def is_valid_cns_file(file: Path) -> bool:
    """Validate that the CNV segment file exists, is non-empty, has a '.cnr' suffix, and has the expected header."""
    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and has_correct_suffix(file, ".cns")
        and has_expected_header(
            file,
            "chromosome\tstart\tend\tgene\tlog2\tdepth\tprobes\tweight\tci_lo\tci_hi",
        )
    )


def is_valid_call_cns_file(file: Path) -> bool:
    """Validate that the call CNV segment file exists, is non-empty, has a '.cnr' suffix, and has the expected header."""
    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and has_correct_suffix(file, ".cns")
        and has_expected_header(
            file,
            "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tp_ttest\tprobes\tweight",
        )
    )


def is_valid_centred_file(file: Path) -> bool:
    """Validate that the call CNV segment file exists, is non-empty, has a '.cnr' suffix, and has the expected header."""
    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and has_correct_suffix(file, "centred.cns")
    )


def is_valid_bintest_cns_file(file: Path) -> bool:
    """Validate that the bintest segmentation file exists, is non-empty, and has the expected suffix."""
    return file_exists_and_nonempty(file) and has_correct_suffix(file, ".bintest.cns")


def is_valid_genemetrics_file(file: Path) -> bool:
    """Validate that the genemetrics file exists, is non-empty, and has an expected suffix."""
    valid_suffixes = (".ratios.genemetrics.out", ".segments.genemetrics.out")
    valid_headers = (
        "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tprobes",
        "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tcn\tp_ttest\tprobes\tsegment_weight\tsegment_probes",
        "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tci_hi\tci_lo\tprobes\tsegment_weight\tsegment_probes",
    )

    return (
        file_exists_and_nonempty(file)
        and is_regular_file(file)
        and any(file.name.endswith(suffix) for suffix in valid_suffixes)
        and any(has_expected_header(file, header) for header in valid_headers)
    )


def is_valid_scatter_plot(file: Path) -> bool:
    """Validate that the scatter plot exists, is non-empty, and has '-scatter.png' as its suffix."""
    return file_exists_and_nonempty(file) and file.name.endswith("-scatter.png")


def is_valid_diagram_plot(file: Path) -> bool:
    """Validate that the diagram plot exists, is non-empty, and has '-diagram.pdf' as its suffix."""
    return file_exists_and_nonempty(file) and file.name.endswith("-diagram.pdf")


def is_valid_filtered_genemetrics_file(file: Path) -> bool:
    """Validate that the filtered genemetrics file exists, is non-empty, and has the expected suffix."""
    return file_exists_and_nonempty(file) and has_correct_suffix(
        file, ".genemetrics.significant.out"
    )


# -----------------------------------------------------------------------------
# CNVKit command wrappers
# -----------------------------------------------------------------------------
def run_cnvkit_access(reference_fasta: Path, outdir: Path) -> Path:
    """Run cnvkit.py access on a given reference FASTA file

    Args:
        reference_fasta (Path): Path to the reference FASTA file
        outdir (Path): Path to directory where output file will be stored

    Returns:
        Path: Path to the output BED file containing sequence-accessible coordinates
    """
    logger.info("Preparing to run cnvkit.py access...")
    if is_fasta(reference_fasta):
        fasta_stem = reference_fasta.stem
        output_bed_name = f"access-{fasta_stem}.bed"
        output_bed_path = outdir / output_bed_name

        # Skip generation if the file already exists and is valid
        if skip_file_generation(output_bed_path, validator=is_bed):
            return output_bed_path

        cmd = f"cnvkit.py access {str(reference_fasta)} -o {str(output_bed_path)}"
        run_command(cmd)
    else:
        raise ValueError(
            f"{str(reference_fasta)} is not a valid FASTA file. Check input data."
        )

    return output_bed_path


def run_cnvkit_autobin(
    bam_files: t.List[Path],
    baitset_bed: Path,
    access_bed: Path,
    refflat_file: Path,
    outdir: Path,
) -> t.Dict[str, Path]:
    """Run cnvkit.py autobin to produce target and antitarget BED files

    Args:
        bam_files (list[Path]): Paths to BAM files
        baitset_bed (Path): Path to the baitset BED file
        access_bed (Path): Path to the BED file generated by cnvkit.py access
        refflat_file (Path): Path to the refFlat gene annotation file
        outdir (Path): Path to directory where output files will be stored

    Returns:
        dict[str, Path]: A dictionary containing the target and antitarget BED files.
                         Keys are 'target' and 'antitarget'.
    """
    logger.info("Preparing to run cnvkit.py autobin...")
    # Validate input baitset BED and construct output BED paths
    if is_bed(baitset_bed):
        output_target_bed_name = baitset_bed.name.replace(".bed", ".target.bed")
        output_antitarget_bed_name = baitset_bed.name.replace(".bed", ".antitarget.bed")
        output_target_bed_path = outdir / output_target_bed_name
        output_antitarget_bed_path = outdir / output_antitarget_bed_name
    else:
        raise ValueError(
            f"{str(baitset_bed)} is not a valid BED file. Check input data"
        )

    # Validate input BAMs and construct BAM file string for cnvkit.py autobin command
    if all(is_bam(file) for file in bam_files):
        bam_files_as_string = map(str, bam_files)
    else:
        raise ValueError("One or more BAM files are invalid.")

    # If both output files already exist and are valid, skip autobin
    if skip_file_generation(
        output_target_bed_path, validator=is_bed
    ) and skip_file_generation(output_antitarget_bed_path, validator=is_bed):
        return {
            "target": output_target_bed_path,
            "antitarget": output_antitarget_bed_path,
        }

    # Construct cnvkit.py autobin command
    cmd = (
        f'cnvkit.py autobin {" ".join(bam_files_as_string)} -t {str(baitset_bed)} '
        f"-g {str(access_bed)} --annotate {str(refflat_file)} "
        f"--target-output-bed {str(output_target_bed_path)} "
        f"--antitarget-output-bed {str(output_antitarget_bed_path)}"
    )

    # Run cnvkit.py autobin command
    run_command(cmd)

    return {"target": output_target_bed_path, "antitarget": output_antitarget_bed_path}


def run_cnvkit_coverage(bam_file: Path, interval_bed: Path, outdir: Path) -> Path:
    """Run cnvkit.py coverage on a BAM file with a given interval BED file.

    Args:
        bam_file (Path): Path to the BAM file.
        interval_bed (Path): Path to the interval BED file.
        outdir (Path): Directory to store the output coverage file.

    Returns:
        Path: The output coverage file path.
    """
    is_valid_flag = True

    # Check if BAM file is valid
    if is_bam(bam_file):
        logger.info(f"BAM file {str(bam_file)} is valid")
    else:
        logger.error(f"BAM file {str(bam_file)} is not valid. Please check input data.")
        is_valid_flag = False

    # Check if BED file is valid
    if is_bed(interval_bed):
        logger.info(f"BED file {str(interval_bed)} is valid")
    else:
        logger.error(
            f"BED file {str(interval_bed)} is not valid. Please check input data."
        )
        is_valid_flag = False

    if not is_valid_flag:
        raise ValueError(
            f"{str(bam_file)} is not a valid BAM file and/or {str(interval_bed)} is not a valid BED file. Please check input data."
        )

    target_type = interval_bed.stem.split(".")[-1]
    output_file_name = bam_file.name.replace(".bam", f".{target_type}coverage.cnn")
    output_file_path = outdir / output_file_name

    # Skip generation if the coverage file exists and is valid
    if skip_file_generation(output_file_path, validator=is_valid_coverage_file):
        return output_file_path

    cmd = f"cnvkit.py coverage {str(bam_file)} {str(interval_bed)} -o {str(output_file_path)}"
    run_command(cmd)

    return output_file_path


def run_cnvkit_reference(
    coverage_files: t.List[Path],
    reference_fasta: Path,
    output_prefix: str,
    outdir: Path,
    sex: t.Literal["male", "female"],
) -> Path:
    """Run cnvkit.py reference to generate a copy number reference.

    Args:
        coverage_files (List[Path]): List of coverage file paths.
        reference_fasta (Path): Path to the reference FASTA.
        output_prefix (str): Prefix for the output file.
        outdir (Path): Directory to store the output file.
        sex (str): 'male' or 'female'.

    Returns:
        Path: Path to the generated reference file.
    """
    is_valid_flag = True

    if is_fasta(reference_fasta):
        logger.info(f"Reference FASTA file {str(reference_fasta)} is valid.")
    else:
        is_valid_flag = False
        raise ValueError("Please provide a valid reference FASTA file.")

    if is_valid_flag:
        coverage_files_as_string = convert_file_list_to_string(coverage_files)
        output_file_name = f"{output_prefix}.reference.cnn"
        output_file_path = outdir / output_file_name

        cnvkit_reference_command = f"cnvkit.py reference {coverage_files_as_string} -f {str(reference_fasta)} -o {str(output_file_path)}"
        if sex == "male":
            cnvkit_reference_command += " -y"

        if skip_file_generation(output_file_path, validator=is_valid_reference_file):
            return output_file_path

        run_command(cnvkit_reference_command)

        return output_file_path


def run_cnvkit_fix(
    target_coverage_file: Path,
    antitarget_coverage_file: Path,
    copy_number_reference_file: Path,
    output_prefix: str,
    outdir: Path,
) -> Path:
    """Run cnvkit.py fix to generate a copy number ratio file.

    Args:
        target_coverage_file (Path): Path to target coverage file.
        antitarget_coverage_file (Path): Path to antitarget coverage file.
        copy_number_reference_file (Path): Path to the copy number reference.
        output_prefix (str): Prefix for the output file.
        outdir (Path): Directory to store the output file.

    Returns:
        Path: Path to the generated CNV ratio file.
    """
    output_file_name = f"{output_prefix}.cnr"
    output_file_path = outdir / output_file_name

    if skip_file_generation(output_file_path, validator=is_valid_cnr_file):
        return output_file_path

    cnvkit_fix_command = (
        f"cnvkit.py fix {str(target_coverage_file)} {str(antitarget_coverage_file)} "
        f"{str(copy_number_reference_file)} -o {str(output_file_path)}"
    )

    run_command(cnvkit_fix_command)

    return output_file_path


def run_cnvkit_genemetrics(
    ratio_file: Path,
    threshold: float,
    min_probes: int,
    output_prefix: str,
    outdir: Path,
    sex: t.Literal["male", "female"],
    segment_file: t.Optional[Path] = None,
) -> Path:
    """
    Run cnvkit.py genemetrics on a given bin-level log2 ratio file.
    Optionally include a segmentation file. Returns the path to the output genemetrics file.

    Parameters:
        ratio_file (Path): Path to the ratio file (Sample.cnr).
        threshold (float): Threshold value for calling copy number variations.
        min_probes (int): Minimum number of probes required for a call.
        output_prefix (str): Prefix for the output file.
        outdir (Path): Directory to save the output file.
        sex (str): 'male' or 'female'.
        segment_file (Path, optional): Path to the segmentation file (Sample.cns).

    Returns:
        Path: Path to the generated genemetrics output file.
    """
    suffix = "segments" if segment_file else "ratios"
    output_file_name = f"{output_prefix}.{suffix}.genemetrics.out"
    output_file_path = outdir / output_file_name

    if skip_file_generation(output_file_path, validator=is_valid_genemetrics_file):
        return output_file_path

    cnvkit_genemetrics_command = (
        f"cnvkit.py genemetrics {ratio_file} -t {threshold} -m {min_probes} "
        f"-o {output_file_path}"
    )

    if segment_file:
        cnvkit_genemetrics_command += f" -s {segment_file}"
    if sex == "male":
        cnvkit_genemetrics_command += " -y"

    run_command(cnvkit_genemetrics_command)

    return output_file_path


def _validate_output_files_for_bam(
    bam: Path,
    outdir: Path,
    expected_files: t.List[t.Tuple[str, t.Callable[[Path], bool]]],
) -> bool:
    """
    Check if all expected cnvkit.py batch output files for a given BAM file exist and are valid.

    Returns:
        bool: True if all files are valid, False otherwise.
    """
    for suffix, validator in expected_files:
        expected_file = outdir / f"{bam.stem}{suffix}"
        if not skip_file_generation(expected_file, validator=validator):
            logger.debug(f"File {expected_file} is missing or invalid.")
            return False
    return True


def run_cnvkit_batch(
    tumour_bams: t.List[Path],
    copy_number_reference_file: Path,
    outdir: Path,
    sex: t.Literal["male", "female"],
) -> Path:
    """
    Run the CNVkit batch command on a list of tumour BAM files.

    This function uses the existing helper functions to check if the expected CNVKit batch
    output files already exist. For each tumour BAM, it expects a corresponding CNV ratio
    file (with the '.cnr' suffix) in the output directory. If all such files exist and are valid,
    the batch run is skipped.

    Args:
        tumour_bams (List[Path]): List of tumour BAM files.
        copy_number_reference_file (Path): Path to the copy number reference.
        outdir (Path): Directory to store the batch results.
        sex (Literal["male", "female"]): Sex of the sample, which affects reference options.

    Returns:
        Path: Path to the output directory containing the CNVKit batch results.
    """
    # Verify R dependencies and set up the environment.
    verify_R_dependencies()
    env = setup_R_environment()

    # Define the expected output file types and validators.
    expected_suffixes_and_validators = [
        (".targetcoverage.cnn", is_valid_targetcoverage_file),
        (".antitargetcoverage.cnn", is_valid_antitargetcoverage_file),
        (".cnr", is_valid_cnr_file),
        (".cns", is_valid_cns_file),
        (".bintest.cns", is_valid_bintest_cns_file),
        (".call.cns", is_valid_call_cns_file),
    ]

    # Check validity of all output files for each tumour BAM.
    outputs_valid = all(
        _validate_output_files_for_bam(bam, outdir, expected_suffixes_and_validators)
        for bam in tumour_bams
    )

    if outputs_valid:
        logger.info(
            "All CNVkit batch output files already exist and are valid. Skipping CNVkit batch run."
        )
        return outdir

    # If outputs are missing or are not valid, continue with the cnvkit.py batch run
    # Build the tumour BAM files string for the batch command.
    tumour_bams_as_string = convert_file_list_to_string(tumour_bams)
    cnvkit_batch_cmd = (
        f"cnvkit.py batch {tumour_bams_as_string} "
        "-m hybrid "
        "--drop-low-coverage "
        f"--reference {str(copy_number_reference_file)} "
        f"--output-dir {str(outdir)}"
    )

    # If the sample sex is male, add the male reference option.
    if sex == "male":
        cnvkit_batch_cmd += " --male-reference"

    # Run the CNVkit batch command.
    run_command(cnvkit_batch_cmd, env=env)

    return outdir


def run_cnvkit_scatter(
    ratio_file: Path, segment_file: Path, output_directory: Path
) -> Path:
    """Run cnvkit.py scatter on a given ratio and segment file. Returns a Path object to the output plot"""
    output_scatter_plot_filename = ratio_file.name.replace(".cnr", "-scatter.png")
    output_scatter_plot_path = output_directory / output_scatter_plot_filename

    cnvkit_scatter_cmd = f"cnvkit.py scatter {str(ratio_file)} -s {str(segment_file)} -o {str(output_scatter_plot_path)}"
    logger.debug(f"Constructed cnvkit.py scatter command: {cnvkit_scatter_cmd}")
    logger.info("Running cnvkit.py scatter command ...")

    if skip_file_generation(output_scatter_plot_path, validator=is_valid_scatter_plot):
        return output_scatter_plot_path

    run_command(cnvkit_scatter_cmd)

    return output_scatter_plot_path


def run_cnvkit_diagram(
    ratio_file: Path, segment_file: Path, output_directory: Path
) -> Path:
    """Run cnvkit.py diagram on a given ratio and segment file. Returns a Path object to the output plot"""
    output_diagram_plot_filename = ratio_file.name.replace(".cnr", "-diagram.pdf")
    output_diagram_plot_path = output_directory / output_diagram_plot_filename

    cnvkit_diagram_cmd = f"cnvkit.py diagram {str(ratio_file)} -s {str(segment_file)} -o {str(output_diagram_plot_path)}"
    logger.debug(f"Constructed cnvkit.py diagram command: {cnvkit_diagram_cmd}")
    logger.info("Running cnvkit.py diagram command ...")

    if skip_file_generation(output_diagram_plot_path, validator=is_valid_diagram_plot):
        return output_diagram_plot_path

    run_command(cnvkit_diagram_cmd)

    return output_diagram_plot_path


# -----------------------------------------------------------------------------
# CNVKit batch postprocessing functions
# -----------------------------------------------------------------------------
def _check_centring_method(centring_method: str):
    # Define the allowed centring methods.
    allowed_methods = {"mean", "median", "mode", "biweight"}
    if centring_method not in allowed_methods:
        raise ValueError(
            f"Invalid centring_method: {centring_method}. "
            f"Must be one of: {', '.join(sorted(allowed_methods))}"
        )


def perform_centring(
    copy_number_call_file: Path,
    outdir: Path,
    centring_method: str = "median",
) -> t.Tuple[Path, t.Optional[float]]:
    """
    Perform centring using CNVKit and extract the log2 shift value with the specified centring method.

    Accepted centring methods are: "mean", "median", "mode", "biweight".

    Args:
        copy_number_call_file (Path): The input file containing CNV calls.
        outdir (Path): The output directory to save the centred file.
        centring_method (str, optional): The centring method to use.
            Must be one of "mean", "median", "mode", or "biweight".
            Defaults to "median".

    Returns:
        Tuple[Path, Optional[float]]: A tuple containing:
            - The path to the centred file.
            - The extracted log2 shift value (None if not found).

    Raises:
        ValueError: If an unsupported centring method is provided.
    """
    # Check that the supplied centring method is a valid method
    _check_centring_method(centring_method)

    logger.info(f"Performing {centring_method} centring on {copy_number_call_file} ...")

    # Construct the output file path using the centring method in the suffix.
    centred_suffix = f".{centring_method}_centred{copy_number_call_file.suffix}"
    centred_file_name = copy_number_call_file.stem + centred_suffix
    centred_file_path = outdir / centred_file_name

    # Construct the centring command with the specified centring method.
    centring_cmd = (
        f"cnvkit.py call -m none {copy_number_call_file} --center {centring_method} "
        f"-o {centred_file_path}"
    )

    # If the output file already exists and is valid, skip generation.
    if skip_file_generation(centred_file_path, validator=is_valid_centred_file):
        return centred_file_path, None

    try:
        result = execute_command(centring_cmd)
        log_success(centring_cmd, result)

        # Attempt to extract the log2 shift value from the command's stderr.
        match = re.search(r"Shifting log2 values by ([\d.-]+)", result.stderr)
        shift_value = float(match.group(1)) if match else None

        if shift_value is not None:
            logger.info(f"Extracted log2 shift value: {shift_value}")
        else:
            logger.warning("No log2 shift value found in the output.")

        return centred_file_path, shift_value

    except subprocess.CalledProcessError as e:
        log_error(centring_cmd, e)
        return centred_file_path, None


def filter_genemetrics_file(
    genemetrics_file: Path, lower_threshold: float, upper_threshold: float, outdir: Path
) -> Path:
    """Filter genemetrics file to only show gene-level log2(FC) calls above/below user-defined thresholds."""

    outdir.mkdir(parents=True, exist_ok=True)

    try:
        df = pd.read_csv(genemetrics_file, sep="\t", dtype={"log2": float})
    except Exception as e:
        raise ValueError(f"Error reading the file: {e}")

    if "log2" not in df.columns:
        raise KeyError("The input file does not contain a 'log2' column.")

    filtered_df = df[(df["log2"] <= lower_threshold) | (df["log2"] >= upper_threshold)]

    filtered_genemetrics_file_name = genemetrics_file.name.replace(
        ".genemetrics.out", ".genemetrics.significant.out"
    )
    filtered_genemetrics_file_path = outdir / filtered_genemetrics_file_name

    if skip_file_generation(
        filtered_genemetrics_file_path, validator=is_valid_filtered_genemetrics_file
    ):
        return filtered_genemetrics_file_path

    filtered_df.to_csv(filtered_genemetrics_file_path, sep="\t", index=False)

    return filtered_genemetrics_file_path


# -----------------------------------------------------------------------------
# CNVKit file parsing and filtering
# -----------------------------------------------------------------------------
def parse_genemetrics_file(file_path: Path):
    """Parse the genemetrics output file and return a list of log2(FC) values"""
    logger.info(f"Parsing genemetrics file {str(file_path)}")
    log2_values = []
    with file_path.open() as f:
        next(f)  # Skip the header line
        for line in f:
            parts = line.strip().split("\t")
            log2_value = float(parts[4])
            log2_values.append(log2_value)

    logger.debug(f"Log2(FC) values from {str(file_path)}: {log2_values}")
    return log2_values


def filter_unplaced_contigs_from_cnvkit_output_file(
    cnvkit_output_file: Path, unplaced_contigs: t.List[str]
):
    logger.info(f"Removing unplaced contigs from {str(cnvkit_output_file)}...")

    with cnvkit_output_file.open() as f:
        lines = f.readlines()

    filtered_lines = [
        line
        for line in lines
        if not any(line.startswith(contig) for contig in unplaced_contigs)
    ]

    with cnvkit_output_file.open("w") as f:
        f.writelines(filtered_lines)

    logger.info(
        f"Unplaced contigs removed from {str(cnvkit_output_file)}. Filtered {len(lines) - len(filtered_lines)} lines. Output saved to {str(cnvkit_output_file)}"
    )


# -----------------------------------------------------------------------------
# R related functions
# -----------------------------------------------------------------------------


def is_R_package_installed(package_name: str) -> bool:
    """Check if an R package is installed.

    Args:
        package_name (str): The name of the R package to check.

    Returns:
        bool: True if the package is installed, False otherwise.
    """
    cmd = f'Rscript -e \'if (!require("{package_name}", quietly = TRUE)) stop("{package_name} package is not installed.")\''
    try:
        execute_command(cmd)
        result = True
    except subprocess.CalledProcessError:
        result = False
    return result


def get_libPaths_used_by_R() -> t.List[str]:
    """Get the library paths used by R.

    Returns:
        List[str]: A list of library paths used by R.
    """
    cmd = "Rscript -e 'cat(.libPaths(), sep = \"\\n\")'"
    result = execute_command(cmd)
    libPaths = result.stdout.strip().split("\n")
    return libPaths


def get_concatenated_libPaths_used_by_R() -> str:
    """Get the concatenated library paths used by R.

    Returns:
        str: A string of concatenated library paths used by R.
    """
    return ":".join(get_libPaths_used_by_R())


def verify_R_dependencies() -> None:
    """Ensure that the necessary R package is installed."""
    if not is_R_package_installed("DNAcopy"):
        msg = "DNAcopy package is not installed. Please install the package."
        logger.error(msg)
        raise ValueError(msg)
    logger.info("DNAcopy package is installed.")


def setup_R_environment() -> t.Dict[str, str]:
    """Set up R environment variables if they are not already set."""
    env = os.environ.copy()
    if not env.get("R_LIBS_USER"):
        r_libs = get_concatenated_libPaths_used_by_R()
        env["R_LIBS_USER"] = r_libs
        logger.info(f"R_LIBS_USER not set. Setting R_LIBS_USER to: {r_libs}")
    return env
