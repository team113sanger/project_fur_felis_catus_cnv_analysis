import pytest
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict
from tempfile import TemporaryDirectory

from tests.mocks.mock_files import get_example_sample_metadata_xlsx
from fur_cnvkit.utils.fur_utils import (
    get_sample_id_from_file_path,
    get_sample_ids_for_file_list,
    get_sample_specific_files,
    validate_metadata_keys,
    map_sample_ids_to_study_ids,
    convert_paths,
    extract_sample_ids_from_exclude_file,
    filter_files_by_exclude_samples,
    remove_unwanted_sample_files,
    determine_sample_sexes,
    split_file_list_by_sample_sex,
    get_tumour_normal_status,
    categorise_files_by_tumour_normal_status,
    get_sample_study,
    set_metadata_columns,
    DEFAULT_METADATA_COLUMNS,
)


# Fixtures


@pytest.fixture
def example_sample_metadata_xlsx():
    sample_metadata_xlsx = get_example_sample_metadata_xlsx()
    return sample_metadata_xlsx


@pytest.fixture(autouse=True)
def reset_metadata_configuration():
    """Ensure metadata column configuration is reset before and after each test."""
    set_metadata_columns(None)
    yield
    set_metadata_columns(None)


# Tests


def test_get_sample_id_from_file_path():
    assert get_sample_id_from_file_path(Path("sample1.bam")) == "sample1"
    assert get_sample_id_from_file_path(Path("sample1.bam.old")) == "sample1"
    assert get_sample_id_from_file_path(Path("sample2.txt")) == "sample2"


def test_get_sample_ids_for_file_list():
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample2.bam"),
        Path("/path/to/sample1.bam"),  # Duplicate
    ]
    expected = {"sample1", "sample2"}
    assert get_sample_ids_for_file_list(files) == expected


def test_get_sample_specific_files_when_suffixes_is_none_and_files_found():
    """
    If suffixes is None, we return all files containing 'sample1' in their filenames.
    """
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample1.txt"),
        Path("/path/to/another_sample1.fastq"),
        Path("/path/to/unrelated_file.txt"),
    ]
    sample = "sample1"

    result = get_sample_specific_files(files, sample, suffixes=None)
    expected = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample1.txt"),
        Path("/path/to/another_sample1.fastq"),
    ]
    assert sorted(result) == sorted(expected)


def test_get_sample_specific_files_when_suffixes_is_none_and_no_files_found():
    """
    If suffixes is None but no file contains the sample string,
    we expect a ValueError.
    """
    files = [
        Path("/path/to/foo.bam"),
        Path("/path/to/bar.txt"),
    ]
    sample = "sample1"

    with pytest.raises(ValueError, match="No files found containing sample"):
        get_sample_specific_files(files, sample, suffixes=None)


def test_get_sample_specific_files_with_single_suffix_all_good():
    """
    If suffixes has a single item, we need at least one file that ends
    with that suffix (and contains the sample name).
    """
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample1.txt"),
        Path("/path/to/sample2.bam"),
    ]
    sample = "sample1"
    suffixes = [".bam"]

    # Expect to find sample1.bam, which ends with .bam and contains "sample1"
    result = get_sample_specific_files(files, sample, suffixes=suffixes)
    expected = [Path("/path/to/sample1.bam")]
    assert result == expected


def test_get_sample_specific_files_with_single_suffix_none_found():
    """
    If we can't find any file for the single suffix, we raise ValueError.
    """
    files = [
        Path("/path/to/sample1.txt"),
        Path("/path/to/sample1.fastq"),
    ]
    sample = "sample1"
    suffixes = [".bam"]

    # None of these files end with .bam, so we expect an error.
    with pytest.raises(
        ValueError,
        match="No files found containing sample 'sample1' with suffix '.bam'",
    ):
        get_sample_specific_files(files, sample, suffixes=suffixes)


def test_get_sample_specific_files_with_multiple_suffixes_success():
    """
    With multiple suffixes, we need at least one file for each suffix
    that contains the sample in its filename.
    """
    files = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample1.txt"),
        Path("/path/to/sample1.fastq"),  # Not in suffixes, but that’s okay
        Path("/path/to/another_sample1.txt"),
    ]
    sample = "sample1"
    suffixes = [".bam", ".txt"]

    # We have at least one .bam (sample1.bam)
    # and at least one .txt (sample1.txt, another_sample1.txt).
    # The function returns the union of those matches.
    result = get_sample_specific_files(files, sample, suffixes=suffixes)
    expected = [
        Path("/path/to/sample1.bam"),
        Path("/path/to/sample1.txt"),
        Path("/path/to/another_sample1.txt"),
    ]
    assert sorted(result) == sorted(expected)


def test_get_sample_specific_files_with_multiple_suffixes_missing_one():
    """
    If any one suffix fails to match any file, we should raise a ValueError immediately.
    """
    files = [
        Path("/path/to/sample1.bam"),  # .bam is there
        Path("/path/to/sample1.fastq"),
        Path("/path/to/sample2.txt"),
    ]
    sample = "sample1"
    suffixes = [".bam", ".txt"]

    # We do have sample1.bam (good for .bam),
    # but there's NO file for sample1 that ends with .txt,
    # so we expect ValueError.
    with pytest.raises(
        ValueError,
        match="No files found containing sample 'sample1' with suffix '.txt'",
    ):
        get_sample_specific_files(files, sample, suffixes=suffixes)


def test_validate_metadata_keys_success(tmp_path):
    """
    Test validate_metadata_keys for a successful scenario when no keys are missing.
    """
    metadata = {
        "all_bams": ["sample1.bam", "sample2.bam"],
        "tumour_bams": ["sample1.bam"],
        "normal_bams": ["sample2.bam"],
        "reference_fasta": "ref.fa",
        "unplaced_contig_prefixes": ["foo"],
        "baitset_bed": "baits.bed",
        "refflat_file": "refflat.txt",
        "sample_metadata_xlsx": "metadata.xlsx",
        "access_bed": "access.bed",
        "targets_bed": "targets.bed",
        "antitargets_bed": "antitargets.bed",
    }
    expected_keys = (
        "all_bams",
        "tumour_bams",
        "normal_bams",
        "reference_fasta",
        "unplaced_contig_prefixes",
        "baitset_bed",
        "refflat_file",
        "sample_metadata_xlsx",
        "access_bed",
        "targets_bed",
        "antitargets_bed",
    )

    # Should not raise an error
    try:
        validate_metadata_keys(metadata, expected_keys, tmp_path / "param.json")
    except KeyError:
        pytest.fail("validate_metadata_keys raised KeyError unexpectedly.")


def test_validate_metadata_keys_missing(tmp_path):
    """
    Test validate_metadata_keys when some keys are missing.
    """
    metadata = {
        "all_bams": ["sample1.bam", "sample2.bam"],
        # missing tumour_bams, normal_bams, etc.
        "reference_fasta": "ref.fa",
    }
    expected_keys = (
        "all_bams",
        "tumour_bams",
        "normal_bams",
        "reference_fasta",
    )
    with pytest.raises(KeyError) as excinfo:
        validate_metadata_keys(metadata, expected_keys, tmp_path / "param.json")

    # Check if the error mentions missing keys
    assert "tumour_bams" in str(excinfo.value)
    assert "normal_bams" in str(excinfo.value)


def test_convert_paths(tmp_path):
    """
    Test that convert_paths correctly converts string paths
    to Path objects for all but the excluded keys.
    """
    metadata = {
        "path_1": str(tmp_path / "file1.txt"),
        "path_2": [str(tmp_path / "file2.txt"), str(tmp_path / "file3.txt")],
        "some_value": 123,
        "exclude_this": str(tmp_path / "exclude_me.txt"),
    }
    exclude_keys = {"exclude_this"}

    converted = convert_paths(metadata, exclude_keys)
    assert isinstance(converted["path_1"], Path)
    assert all(isinstance(p, Path) for p in converted["path_2"])
    # This key should remain unchanged because it's excluded
    assert isinstance(converted["exclude_this"], str)
    # Non-path items should remain as they are
    assert converted["some_value"] == 123


def test_map_sample_ids_to_study_ids_success():
    """
    Test map_sample_ids_to_study_ids with valid sample IDs in the metadata.
    """
    # Create a mock DataFrame with two sheets
    with TemporaryDirectory() as temp_dir:
        metadata_file = Path(temp_dir) / "metadata.xlsx"
        with pd.ExcelWriter(metadata_file) as writer:
            df_study1 = pd.DataFrame(
                {
                    "Sanger DNA ID": ["S1", "S2"],
                    "Some Column": ["val1", "val2"],
                }
            )
            df_study2 = pd.DataFrame(
                {
                    "Sanger DNA ID": ["S3"],
                    "Some Column": ["val3"],
                }
            )
            df_study1.to_excel(writer, sheet_name="Study1", index=False)
            df_study2.to_excel(writer, sheet_name="Study2", index=False)

        sample_ids = {"S1", "S2", "S3"}
        result = map_sample_ids_to_study_ids(sample_ids, metadata_file)

        # We expect each study to have a list of the sample IDs that appear in it
        expected = {
            "Study1": ["S1", "S2"],
            "Study2": ["S3"],
        }
        # Order of items in the lists might differ, so we'll compare sorted lists
        for study in result:
            assert sorted(result[study]) == sorted(expected[study])


def test_map_sample_ids_to_study_ids_with_study_column(tmp_path):
    metadata_df = pd.DataFrame(
        {
            "Sanger_DNA_ID": ["S1", "S2", "S3"],
            "Study": ["StudyA", "StudyB", "StudyA"],
        }
    )
    metadata_file = tmp_path / "metadata.tsv"
    metadata_df.to_csv(metadata_file, sep="\t", index=False)

    set_metadata_columns(
        {
            "sample_id": "Sanger_DNA_ID",
            "tumour_normal": DEFAULT_METADATA_COLUMNS.tumour_normal,
            "sex": DEFAULT_METADATA_COLUMNS.sex,
            "study": "Study",
        }
    )

    result = map_sample_ids_to_study_ids({"S1", "S2", "S3"}, metadata_file)
    expected = {"StudyA": {"S1", "S3"}, "StudyB": {"S2"}}
    assert {study: set(samples) for study, samples in result.items()} == expected


def test_map_sample_ids_to_study_ids_missing_sample():
    """
    Test map_sample_ids_to_study_ids with a sample ID that does not appear
    in any sheet, which should raise ValueError from get_sample_study.
    """
    with TemporaryDirectory() as temp_dir:
        metadata_file = Path(temp_dir) / "metadata.xlsx"
        with pd.ExcelWriter(metadata_file) as writer:
            df_study1 = pd.DataFrame({"Sanger DNA ID": ["S1"], "Col": ["val1"]})
            df_study1.to_excel(writer, sheet_name="Study1", index=False)

        # "S2" won't be found
        sample_ids = {"S1", "S2"}

        with pytest.raises(ValueError, match="Sample ID 'S2' not found in metadata."):
            map_sample_ids_to_study_ids(sample_ids, metadata_file)


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
    expected = {"sample1": "male", "sample2": "female", "sample3": "unknown"}
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

    with pytest.raises(
        ValueError,
        match=f"Sample ID 'sample4' is missing from the metadata file '{excel_file}'.",
    ) as excinfo:
        determine_sample_sexes(files, excel_file)

    expected_error_message = (
        f"Sample ID 'sample4' is missing from the metadata file '{excel_file}'."
    )
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
            "male": [Path("/path/to/sample1.bam"), Path("/path/to/sample4.bam")],
            "female": [Path("/path/to/sample2.bam"), Path("/path/to/sample3.bam")],
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

    with pytest.raises(
        ValueError,
        match=f"Sample ID 'sample3' is missing from the metadata file '{excel_file}'.",
    ) as excinfo:
        split_file_list_by_sample_sex(files, excel_file)

    expected_error_message = (
        f"Sample ID 'sample3' is missing from the metadata file '{excel_file}'."
    )
    assert expected_error_message in str(excinfo.value)
    assert expected_error_message in caplog.text


def test_get_tumour_normal_status_valid_sample_id(example_sample_metadata_xlsx):
    # Test with a valid sample ID (replace 'VALID_SAMPLE_ID' with an actual ID from the file)
    sample_id = "CATD0161a"
    status = get_tumour_normal_status(example_sample_metadata_xlsx, sample_id)
    assert status in ["T", "N"], "The status should be 'T' or 'N'."
    assert status == "T"


def test_get_tumour_normal_status_invalid_sample_id(example_sample_metadata_xlsx):
    # Test with an invalid sample ID
    sample_id = "INVALID_SAMPLE_ID"
    with pytest.raises(
        ValueError, match=f"Sample ID '{sample_id}' not found in metadata."
    ):
        get_tumour_normal_status(example_sample_metadata_xlsx, sample_id)


def test_get_tumour_normal_status_missing_columns():
    # Create a temporary DataFrame with missing columns and save to an Excel file
    df = pd.DataFrame(
        {"Some Other Column": ["A", "B", "C"], "Another Column": [1, 2, 3]}
    )
    temp_file = Path("/tmp/missing_columns.xlsx")
    df.to_excel(temp_file, index=False, sheet_name="Sheet1")

    with pytest.raises(ValueError) as excinfo:
        get_tumour_normal_status(temp_file, "ANY_ID")
    assert "Required column(s) Sanger DNA ID, T/N are missing" in str(excinfo.value)

    # Clean up the temporary file
    temp_file.unlink()


def test_get_tumour_normal_status_multiple_patient_samples():
    # Create a DataFrame to test single character extraction
    df = pd.DataFrame(
        {
            "Sanger DNA ID": ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"],
            "T/N": ["T1", "N2", "T", "N1", "T5"],  # Includes different formats of T/N
        }
    )
    temp_file = Path("/tmp/test_single_character_code.xlsx")
    with pd.ExcelWriter(temp_file) as writer:
        df.to_excel(writer, sheet_name="Sheet1", index=False)

    try:
        # Assert the single character extraction
        assert get_tumour_normal_status(temp_file, "Sample1") == "T"
        assert get_tumour_normal_status(temp_file, "Sample2") == "N"
        assert get_tumour_normal_status(temp_file, "Sample3") == "T"
        assert get_tumour_normal_status(temp_file, "Sample4") == "N"
        assert get_tumour_normal_status(temp_file, "Sample5") == "T"

    finally:
        # Cleanup
        temp_file.unlink()


def test_get_tumour_normal_status_with_custom_columns(tmp_path):
    df = pd.DataFrame(
        {
            "Sanger_DNA_ID": ["Sample1", "Sample2"],
            "Phenotype": ["Tumour", "Normal"],
            "Sex_Column": ["M", "F"],
        }
    )
    tsv_file = tmp_path / "metadata.tsv"
    df.to_csv(tsv_file, sep="\t", index=False)

    set_metadata_columns(
        {
            "sample_id": "Sanger_DNA_ID",
            "tumour_normal": "Phenotype",
            "sex": "Sex_Column",
        }
    )

    assert get_tumour_normal_status(tsv_file, "Sample1") == "T"
    assert get_tumour_normal_status(tsv_file, "Sample2") == "N"
    assert determine_sample_sexes(
        [Path("Sample1.bam"), Path("Sample2.bam")], tsv_file
    ) == {"Sample1": "male", "Sample2": "female"}


def test_get_tumour_normal_status_invalid_status():
    # Create a DataFrame to test single character extraction, including an invalid T/N status
    df = pd.DataFrame(
        {
            "Sanger DNA ID": [
                "Sample1",
                "Sample2",
                "Sample3",
                "Sample4",
                "Sample5",
                "SampleInvalid",
            ],
            "T/N": [
                "T1",
                "N2",
                "T",
                "N1",
                "T5",
                "InvalidStatus",
            ],  # Includes valid and invalid T/N statuses
        }
    )
    temp_file = Path("/tmp/test_single_character_code_with_invalid.xlsx")
    with pd.ExcelWriter(temp_file) as writer:
        df.to_excel(writer, sheet_name="Sheet1", index=False)

    try:
        # Test invalid T/N status
        with pytest.raises(
            ValueError,
            match="Invalid T/N status 'InvalidStatus' for sample ID 'SampleInvalid'.",
        ):
            get_tumour_normal_status(temp_file, "SampleInvalid")

    finally:
        # Cleanup
        temp_file.unlink()


def test_categorise_files_by_tumour_normal_status():
    # Mock sample metadata
    df = pd.DataFrame(
        {
            "Sanger DNA ID": ["Sample1", "Sample2", "Sample3", "Sample4"],
            "T/N": ["T1", "N2", "T", "N1"],
        }
    )
    temp_metadata_file = Path("/tmp/sample_metadata.xlsx")
    with pd.ExcelWriter(temp_metadata_file) as writer:
        df.to_excel(writer, sheet_name="Sheet1", index=False)

    # Mock files with sample IDs embedded
    files = [
        Path("/path/to/Sample1.fastq"),
        Path("/path/to/Sample2.fastq"),
        Path("/path/to/Sample3.fastq"),
        Path("/path/to/Sample4.fastq"),
    ]

    # Call the function
    result = categorise_files_by_tumour_normal_status(files, temp_metadata_file)

    # Expected output
    expected = defaultdict(list)
    expected["T"] = [Path("/path/to/Sample1.fastq"), Path("/path/to/Sample3.fastq")]
    expected["N"] = [Path("/path/to/Sample2.fastq"), Path("/path/to/Sample4.fastq")]

    # Assert the result matches the expected output
    assert result == expected

    # Cleanup
    temp_metadata_file.unlink()


def test_get_sample_study():
    # Create a temporary directory for the test Excel file
    with TemporaryDirectory() as temp_dir:
        temp_file = Path(temp_dir) / "sample_metadata.xlsx"

        # Create a sample Excel file with multiple sheets
        with pd.ExcelWriter(temp_file) as writer:
            df1 = pd.DataFrame(
                {"Sanger DNA ID": ["S1", "S2"], "Other Info": ["Info1", "Info2"]}
            )
            df2 = pd.DataFrame(
                {"Sanger DNA ID": ["S3", "S4"], "Other Info": ["Info3", "Info4"]}
            )

            df1.to_excel(writer, index=False, sheet_name="Study1")
            df2.to_excel(writer, index=False, sheet_name="Study2")

        # Test cases
        assert get_sample_study(temp_file, "S1") == "Study1"
        assert get_sample_study(temp_file, "S3") == "Study2"

        # Test case for a missing sample ID
        with pytest.raises(ValueError, match="Sample ID 'S5' not found in metadata."):
            get_sample_study(temp_file, "S5")
