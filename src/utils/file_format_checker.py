import logging
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


def is_fasta(file_path: Path) -> bool:
    """Checks if a file is in FASTA format by validating its suffix and content.

    Args:
        file_path (Path): Path to the file.

    Returns:
        bool: True if the file is a valid FASTA file, False otherwise.
    """
    if not file_path.suffix.lower() in {".fasta", ".fa"}:
        logger.warning(f"File {file_path} does not have a valid FASTA suffix.")
        return False

    try:
        with file_path.open("r") as file:
            return any(line.startswith(">") for line in file)
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return False


def is_bam(file_path: Path) -> bool:
    """Validates if the given file is in BAM format. Uses pysam to check BAM validity.

    Args:
        file_path (Path): Path to file to be checked.

    Returns:
        bool:
            - True if the file has a valid BAM suffix (.bam) and can be read by pysam.
            - False if the file does not have a valid BAM suffix, or cannot be read by pysam.
    """
    if not str(file_path).lower().endswith(".bam"):
        logger.warning(f"File {file_path} does not have a valid BAM suffix (.bam).")
        return False

    try:
        with pysam.AlignmentFile(file_path, "rb") as bam_file:
            # Attempt to read the first record to verify BAM format
            for _ in bam_file:
                return True
    except Exception as e:
        logger.error(f"Error validating BAM file {file_path}: {e}")
    return False


def _is_bed(file_path: Path) -> bool:
    """Check if a file has a valid BED format.

    Args:
        file_path (Path): Path to the file to validate.

    Returns:
        bool: True if the file appears to follow the BED format, False otherwise.
    """
    try:
        with open(file_path, "r") as file:
            lines = [
                line.strip()
                for line in file
                if line.strip() and not line.strip().startswith("#")
            ]  # Extract all non-empty and non-comment lines

        for line in lines:
            fields = line.split("\t")
            expected_minimum_fields = 3
            is_start_column_digit = fields[1].isdigit()
            is_end_column_digit = fields[2].isdigit()
            if (
                len(fields) < expected_minimum_fields
                or not is_start_column_digit
                or not is_end_column_digit
            ):
                return False
        return True
    except Exception as e:
        print(f"Error reading file: {e}")
        return False


def is_bed(file_path: Path) -> bool:
    """Validates if the given file is in BED format. Checks if each line has at least three tab-separated fields.

    Args:
        file_path (Path): Path to file to be checked

    Returns:
        bool:
            - True if file has valid BED suffix (.bed) and has at least three tab-seperated fields
            - False if the file does not have a valid BAM suffix, does not have at least three tab-seperated fields, or cannot be read
    """
    if not str(file_path).lower().endswith(".bed"):
        logger.warning(f"File {file_path} does not have a valid BED suffix (.bed).")
        return False

    try:
        return _is_bed(file_path)
    except Exception as e:
        logger.error(f"Error reading BED file {file_path}: {e}")
        return False


def validate_bam_files(files: list[Path]) -> list[Path]:
    """Validate if a given list of files are in BAM format. Raise an exception is any file is not a valid BAM. Returns a list of paths to valid BAM files.

    Args:
        files (list[Path]): Paths to files to be validated

    Raises:
        ValueError: A given file is not in valid BAM format

    Returns:
        list[Path]: Paths to the validated BAM files
    """
    logger.info("Checking files for valid BAM format...")
    if len(files) == 0:
        raise ValueError("Recieved no files to validate. Check input data.")
    else:
        logger.debug(f"Recieved {len(files)} to validate: {files}")

    valid_bams = []

    for file in files:
        if is_bam(file):
            valid_bams.append(file)
        else:
            raise ValueError(f"{str(file)} is not a valid BAM file. Check input data.")

    logger.info("Successfully validated BAM files")
    logger.debug(f"Valid BAMs: {valid_bams}")

    return valid_bams
