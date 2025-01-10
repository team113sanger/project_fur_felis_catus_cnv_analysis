import pytest
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict

# Import the functions from your script
from utils.fur_utils import (
    get_sample_id_from_file_path,
    get_sample_ids_for_file_list,
    extract_metadata_files_from_config_json,
    extract_sample_ids_from_exclude_file,
    filter_files_by_exclude_samples,
    remove_unwanted_sample_files,
    determine_sample_sexes,
    split_file_list_by_sample_sex,
)


def test_get_sample_id_from_file_path():
    file = Path("/path/to/sample1.bam")
    assert get_sample_id_from_file_path(file) == "sample1"

    file = Path("/another/path/sample2.fastq.gz")
    assert get_sample_id_from_file_path(file) == "sample2"

    file = Path("sample3")
    assert get_sample_id_from_file_path(file) == "sample3"


def test_get_sample_ids_for_file_list():
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample1.bam"),  # Duplicate
    ]
    expected = {"sample1", "sample2"}
    assert get_sample_ids_for_file_list(files) == expected


def test_extract_metadata_files_from_config_json(tmp_path):
    config = {
        "reference_fasta": "/path/to/reference.fasta",
        "baitset_bed": "/path/to/baitset.bed",
        "refflat_file": "/path/to/refflat.txt",
        "sample_metadata_xlsx": "/path/to/metadata.xlsx",
        "targets_bed": "/path/to/targets.bed",
        "antitargets_bed": "/path/to/antitargets.bed",
    }

    config_file = tmp_path / "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f)

    metadata = extract_metadata_files_from_config_json(config_file)
    assert metadata == config


def test_extract_metadata_files_from_config_json_missing_keys(tmp_path):
    config = {
        "reference_fasta": "/path/to/reference.fasta",
        # Missing "baitset_bed" and others
    }

    config_file = tmp_path / "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f)

    with pytest.raises(KeyError) as excinfo:
        extract_metadata_files_from_config_json(config_file)

    assert "baitset_bed" in str(excinfo.value)


def test_extract_sample_ids_from_exclude_file(tmp_path):
    exclude_content = "sample1\nsample3\nsample5"
    exclude_file = tmp_path / "exclude.txt"
    exclude_file.write_text(exclude_content)

    exclude_samples = extract_sample_ids_from_exclude_file(exclude_file)
    assert exclude_samples == ["sample1", "sample3", "sample5"]


def test_extract_sample_ids_from_exclude_file_not_found(tmp_path, caplog):
    exclude_file = tmp_path / "nonexistent.txt"

    with pytest.raises(FileNotFoundError):
        extract_sample_ids_from_exclude_file(exclude_file)

    assert f"Exclude file is not found: {exclude_file}" in caplog.text


def test_filter_files_by_exclude_samples():
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),
        Path("/path/to/another_sample1.bam"),
    ]
    exclude_samples = ["sample1", "sample4"]

    filtered = filter_files_by_exclude_samples(files, exclude_samples)
    expected = [
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),
    ]
    assert filtered == expected


def test_remove_unwanted_sample_files(tmp_path):
    # Create sample files
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),
    ]

    # Create exclude file
    exclude_content = "sample1\nsample3"
    exclude_file = tmp_path / "exclude.txt"
    exclude_file.write_text(exclude_content)

    # Mock the Path objects (assuming actual files are not needed for this test)
    filtered = remove_unwanted_sample_files(files, exclude_file)
    expected = [Path("/path/to/sample2.bam")]
    assert filtered == expected


def test_remove_unwanted_sample_files_no_exclusions(tmp_path):
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
    ]

    exclude_content = ""
    exclude_file = tmp_path / "exclude.txt"
    exclude_file.write_text(exclude_content)

    filtered = remove_unwanted_sample_files(files, exclude_file)
    expected = files.copy()
    assert filtered == expected


def test_determine_sample_sexes(tmp_path):
    # Create sample files
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),
    ]

    # Create a sample metadata Excel file
    metadata = {
        "Study1": {"Sanger DNA ID": ["sample1", "sample2"], "Sex": ["M", "F"]},
        "Study2": {"Sanger DNA ID": ["sample3"], "Sex": ["U"]},
    }
    excel_file = tmp_path / "metadata.xlsx"
    with pd.ExcelWriter(excel_file) as writer:
        for sheet, data in metadata.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sheet, index=False)

    sample_sex_dict = determine_sample_sexes(files, excel_file)
    expected = {"sample1": "M", "sample2": "F", "sample3": "U"}
    assert sample_sex_dict == expected


def test_determine_sample_sexes_invalid_sex(tmp_path, caplog):
    # Create sample files
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
    ]

    # Create a sample metadata Excel file with an invalid sex value
    metadata = {
        "Study1": {
            "Sanger DNA ID": ["sample1", "sample2"],
            "Sex": ["M", "Unknown"],  # "Unknown" is invalid
        }
    }
    excel_file = tmp_path / "metadata.xlsx"
    with pd.ExcelWriter(excel_file) as writer:
        for sheet, data in metadata.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sheet, index=False)

    with pytest.raises(ValueError) as excinfo:
        determine_sample_sexes(files, excel_file)

    assert "Unexpected sex value 'Unknown' for sample 'sample2'" in str(excinfo.value)
    assert "Unexpected sex value 'Unknown' for sample 'sample2'" in caplog.text


def test_determine_sample_sexes_missing_samples(tmp_path, caplog):
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample4.bam"),  # Not in metadata
    ]

    metadata = {"Study1": {"Sanger DNA ID": ["sample1", "sample2"], "Sex": ["M", "F"]}}
    excel_file = tmp_path / "metadata.xlsx"
    with pd.ExcelWriter(excel_file) as writer:
        for sheet, data in metadata.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sheet, index=False)

    with pytest.raises(ValueError) as excinfo:
        determine_sample_sexes(files, excel_file)

    expected_error_message = f"The following samples are missing from the metadata Excel '{excel_file}': sample4."
    assert expected_error_message in str(excinfo.value)
    assert expected_error_message in caplog.text


def test_split_file_list_by_sample_sex(tmp_path):
    # Create sample files
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),
        Path("/path/to/sample4.bam"),
    ]

    # Create a sample metadata Excel file
    metadata = {
        "Study1": {"Sanger DNA ID": ["sample1", "sample2"], "Sex": ["M", "F"]},
        "Study2": {"Sanger DNA ID": ["sample3", "sample4"], "Sex": ["F", "M"]},
    }
    excel_file = tmp_path / "metadata.xlsx"
    with pd.ExcelWriter(excel_file) as writer:
        for sheet, data in metadata.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sheet, index=False)

    file_sex_dict = split_file_list_by_sample_sex(files, excel_file)

    expected = defaultdict(
        list,
        {
            "M": [Path("/path/to/sample1.bam"), Path("/path/to/sample4.bam")],
            "F": [Path("/path/to/sample2.bam"), Path("/path/to/sample3.bam")],
        },
    )

    assert file_sex_dict == expected


def test_split_file_list_by_sample_sex_with_missing_sex(tmp_path, caplog):
    # Create sample files
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample3.bam"),  # Missing sex
    ]

    # Create a sample metadata Excel file
    metadata = {
        "Study1": {"Sanger DNA ID": ["sample1", "sample2"], "Sex": ["M", "F"]}
        # sample3 is missing
    }
    excel_file = tmp_path / "metadata.xlsx"
    with pd.ExcelWriter(excel_file) as writer:
        for sheet, data in metadata.items():
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name=sheet, index=False)

    with pytest.raises(ValueError) as excinfo:
        split_file_list_by_sample_sex(files, excel_file)

    expected_error_message = f"The following samples are missing from the metadata Excel '{excel_file}': sample3."
    assert expected_error_message in str(excinfo.value)
    assert expected_error_message in caplog.text
