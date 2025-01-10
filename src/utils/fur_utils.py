from collections import defaultdict
import json
import logging
from pathlib import Path
import typing as t

import pandas as pd

logger = logging.getLogger(__name__)

# -----------------------------------------------
# Functions for handling sample IDs
# -----------------------------------------------


def get_sample_id_from_file_path(file: Path) -> str:
    sample_id = file.name.split(".")[0]

    return sample_id


def get_sample_ids_for_file_list(file_list: t.List[Path]) -> set:
    sample_ids = set()

    for file in file_list:
        sample_id = get_sample_id_from_file_path(file)
        sample_ids.add(sample_id)

    return sample_ids


# -----------------------------------------------
# Functions for processing the config file
# -----------------------------------------------


def extract_metadata_files_from_config_json(config_json: Path) -> t.Tuple[Path]:
    """Extract all metadata file paths from the config JSON and return as a tuple"""

    logging.info(f"Extracting metadata file paths from {str(config_json)} ...")
    with open(config_json, "r") as json_file:
        metadata = json.load(json_file)

    # Define the expected keys from the config JSON
    expected_keys = (
        "reference_fasta",
        "baitset_bed",
        "refflat_file",
        "sample_metadata_xlsx",
        "targets_bed",
        "antitargets_bed",
    )

    # Check that the expected keys are present
    missing_keys = tuple(key for key in expected_keys if key not in metadata)

    if missing_keys:
        raise KeyError(
            f"The following keys are missing in {str(config_json)}: {', '.join(missing_keys)}. Please check input data."
        )

    return metadata


# -----------------------------------------------
# Functions for filtering based on the exclude file
# -----------------------------------------------


def extract_sample_ids_from_exclude_file(exclude_file: Path) -> t.List[str]:
    """Read in sample IDs from the exclude file and return as a list"""
    try:
        with open(exclude_file, "r") as file:
            exclude_file_lines = file.readlines()

        exclude_samples = [line.strip() for line in exclude_file_lines]

        return exclude_samples
    except FileNotFoundError:
        logging.error(
            f"Exclude file is not found: {exclude_file}. Please provide an exclude file."
        )
        raise


def filter_files_by_exclude_samples(
    files: t.List[Path], exclude_samples: t.List[str]
) -> t.List[Path]:
    """Remove any files belonging to an exclude sample from the final file list"""
    filtered_files = [
        file
        for file in files
        if not any(sample in file.name for sample in exclude_samples)
    ]

    return filtered_files


def remove_unwanted_sample_files(
    files: t.List[Path], exclude_file: Path
) -> t.List[Path]:
    """Remove files that belong to samples in the exclude file from the final BAM list"""

    # Read in sample IDs from the exlude file
    logging.info(f"Extracting samples to exclude from {str(exclude_file)} ...")
    exclude_samples = extract_sample_ids_from_exclude_file(exclude_file)

    logging.debug(f"Exclude samples: {exclude_samples}")

    # Filter out any files that belong to unwanted samples
    logging.info("Filtering out samples found in the exclude file ...")
    filtered_files = filter_files_by_exclude_samples(files, exclude_samples)

    logging.info(
        f"Successfully filtered out {len(files) - len(filtered_files)} files from unwanted samples, leaving a total of {len(filtered_files)} files."
    )
    logging.debug(f"Filtered files: {filtered_files}")

    return filtered_files


# -----------------------------------------------
# Functions for filtering by sample sex
# -----------------------------------------------


def determine_sample_sexes(
    file_list: t.List[Path], sample_metadata_xlsx: Path
) -> t.Dict[Path, str]:
    # Get the corresponding sample IDs for each file in the file list
    sample_ids = get_sample_ids_for_file_list(file_list)

    # Initialise a dictionary to store each sample's sex
    sample_sex_dict = dict()

    # Load in the sample metadata spreadsheet
    sample_metadata_spreadsheet = pd.ExcelFile(sample_metadata_xlsx)

    # Iterate through each sheet (study) in the metadata spreadsheet
    for sheet_name in sample_metadata_spreadsheet.sheet_names:
        sheet_data = pd.read_excel(sample_metadata_spreadsheet, sheet_name=sheet_name)

        # Filter for rows that include samples from the set of sample ID
        filtered_data = sheet_data[sheet_data["Sanger DNA ID"].isin(sample_ids)]

        # Add the sample name and it's sex to the dictionary
        for _, row in filtered_data.iterrows():
            sample_sex_dict[row["Sanger DNA ID"]] = row["Sex"]

    return sample_sex_dict


def split_file_list_by_sample_sex(
    file_list: t.List[Path], sample_metadata_xlsx: Path
) -> t.DefaultDict[str, t.List[Path]]:
    # Determine the sexes of all samples in the file list
    sample_sex_dict = determine_sample_sexes(file_list, sample_metadata_xlsx)

    # Initialise a dictionary to store sex-seperated files
    file_sex_dict = defaultdict(list)

    for file in file_list:
        sample_id = get_sample_id_from_file_path(file)
        sex = sample_sex_dict.get(sample_id)
        file_sex_dict[sex].append(file)

    return file_sex_dict
