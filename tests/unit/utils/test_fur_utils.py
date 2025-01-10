import pytest
from pathlib import Path
import logging
from unittest.mock import patch, mock_open
from textwrap import dedent

# Import the functions to be tested
from utils.fur_utils import (
    extract_sample_ids_from_exclude_file,
    filter_files_by_exclude_samples,
    remove_unwanted_sample_files,
)


def test_extract_sample_ids_from_exclude_file():
    exclude_file_content = dedent(
        """\
        sample1
        sample2
        sample3
    """
    )

    with patch(
        "builtins.open", mock_open(read_data=exclude_file_content)
    ) as mocked_file:
        result = extract_sample_ids_from_exclude_file(Path("exclude.txt"))

    assert result == ["sample1", "sample2", "sample3"]
    mocked_file.assert_called_once_with(Path("exclude.txt"), "r")


def test_extract_sample_ids_from_exclude_file_file_not_found():
    with pytest.raises(FileNotFoundError):
        extract_sample_ids_from_exclude_file(Path("nonexistent.txt"))


def test_filter_files_by_exclude_samples():
    files = [
        Path("sample1_file1.bam"),
        Path("sample2_file2.bam"),
        Path("sample3_file3.bam"),
        Path("sample4_file4.bam"),
    ]
    exclude_samples = ["sample1", "sample3"]

    result = filter_files_by_exclude_samples(files, exclude_samples)

    assert result == [Path("sample2_file2.bam"), Path("sample4_file4.bam")]


def test_remove_unwanted_sample_files():
    files = [
        Path("sample1_file1.bam"),
        Path("sample2_file2.bam"),
        Path("sample3_file3.bam"),
        Path("sample4_file4.bam"),
    ]
    exclude_file_content = dedent(
        """\
        sample1
        sample3
    """
    )

    with patch(
        "builtins.open", mock_open(read_data=exclude_file_content)
    ) as mocked_file:
        with patch("logging.info") as mock_logging_info:
            result = remove_unwanted_sample_files(files, Path("exclude.txt"))

    assert result == [Path("sample2_file2.bam"), Path("sample4_file4.bam")]
    mocked_file.assert_called_once_with(Path("exclude.txt"), "r")

    mock_logging_info.assert_any_call(
        "Extracting samples to exclude from exclude.txt ..."
    )
    mock_logging_info.assert_any_call(
        "Filtering out samples found in the exclude file ..."
    )
    mock_logging_info.assert_any_call(
        "Successfully filtered out 2 files from unwanted samples, leaving a total of 2 files."
    )
