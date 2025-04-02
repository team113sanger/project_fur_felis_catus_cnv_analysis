import pytest
from pathlib import Path
from unittest.mock import mock_open, patch
from fur_cnvkit.utils.file_format_checker import (
    is_fasta,
    is_bam,
    is_bed,
    validate_bam_files,
)


# Tests for is_fasta
def test_is_fasta_valid(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".fasta")
    content = ">header\nATCG\n"
    mock_file_path.write_text(content)

    assert is_fasta(mock_file_path) is True


def test_is_fasta_no_valid_suffix(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".txt")
    content = ">header\nATCG\n"
    mock_file_path.write_text(content)

    assert is_fasta(mock_file_path) is False


def test_is_fasta_no_header(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".fasta")
    content = "ATCG\n"
    mock_file_path.write_text(content)

    assert is_fasta(mock_file_path) is False


def test_is_fasta_error_reading_file(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".fasta")
    with patch("builtins.open", side_effect=Exception("Error")):
        assert is_fasta(mock_file_path) is False


# Tests for is_bam
@patch("pysam.AlignmentFile")
def test_is_bam_valid(mock_alignment_file, mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".bam")
    mock_alignment_file.return_value.__enter__.return_value = iter(["record"])
    assert is_bam(mock_file_path) is True


@patch("pysam.AlignmentFile")
def test_is_bam_no_valid_suffix(mock_alignment_file, mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".txt")
    assert is_bam(mock_file_path) is False


@patch("pysam.AlignmentFile")
def test_is_bam_invalid_file(mock_alignment_file, mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".bam")
    mock_alignment_file.side_effect = Exception("Error")
    assert is_bam(mock_file_path) is False


# Tests for is_bed
def test_is_bed_valid(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".bed")
    content = "chr1\t100\t200\nchr2\t150\t250\n"
    with patch("builtins.open", mock_open(read_data=content)):
        assert is_bed(mock_file_path) is True


def test_is_bed_no_valid_suffix(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".txt")
    content = "chr1\t100\t200\n"
    with patch("builtins.open", mock_open(read_data=content)):
        assert is_bed(mock_file_path) is False


def test_is_bed_invalid_format(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".bed")
    content = "chr1\t100\n"
    with patch("builtins.open", mock_open(read_data=content)):
        assert is_bed(mock_file_path) is False


def test_is_bed_error_reading_file(mock_file_path):
    mock_file_path = mock_file_path.with_suffix(".bed")
    with patch("builtins.open", side_effect=Exception("Error")):
        assert is_bed(mock_file_path) is False


# Tests for validate_bam_files
@pytest.fixture
def mock_is_bam():
    with patch("fur_cnvkit.utils.file_format_checker.is_bam") as mock_is_bam:
        yield mock_is_bam


def test_validate_bam_files_all_valid(mock_is_bam):
    # Mock is_bam to return True for all files
    mock_is_bam.side_effect = lambda file: True

    files = [Path("file1.bam"), Path("file2.bam")]
    valid_files = validate_bam_files(files)

    assert valid_files == files
    assert mock_is_bam.call_count == len(files)


def test_validate_bam_files_some_invalid(mock_is_bam):
    # Mock is_bam to return False for some files
    mock_is_bam.side_effect = lambda file: file.name == "file1.bam"

    files = [Path("file1.bam"), Path("file2.bam")]

    with pytest.raises(ValueError, match="file2.bam is not a valid BAM file"):
        validate_bam_files(files)

    assert mock_is_bam.call_count == len(files)


def test_validate_bam_files_no_valid_files(mock_is_bam):
    # Mock is_bam to return False for all files
    mock_is_bam.side_effect = lambda file: False

    files = [Path("file1.bam"), Path("file2.bam")]

    with pytest.raises(ValueError, match="file1.bam is not a valid BAM file"):
        validate_bam_files(files)

    # Expecting is_bam to be called only once since the function raises ValueError after the first invalid file
    assert mock_is_bam.call_count == 1


def test_validate_bam_files_empty_list(mock_is_bam):
    # No files to validate
    files = []

    with pytest.raises(
        ValueError, match="Recieved no files to validate. Check input data."
    ):
        validate_bam_files(files)
