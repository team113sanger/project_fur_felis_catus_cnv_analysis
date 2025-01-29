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


def get_sample_specific_files(
    files: t.List[Path], sample: str, suffixes: t.Optional[t.List[str]] = None
) -> t.List[Path]:
    """
    Get the files corresponding to a particular sample.
    - If `suffixes` is None, return all files containing `sample` in their filename.
    - If `suffixes` is a list, then for each suffix in that list, we must find at
      least one file that contains `sample` in its name and ends with that suffix.
      If any suffix doesn't match at least one file, raise a ValueError.

    Args:
        files (List[Path]): List of file paths to filter.
        sample (str): A string (e.g. sample ID) to look for in the filename.
        suffixes (Optional[List[str]]): If provided, we attempt to find files
            for each suffix in this list.

    Returns:
        List[Path]: A list of *all* matching files (for all suffixes) combined.

    Raises:
        ValueError: If no files match the sample (when suffixes is None),
            or if any suffix in `suffixes` has no matching file.
    """
    suffixes_as_string = ", ".join(suffixes if suffixes is not None else [])
    file_message = (
        f"Finding {suffixes_as_string} files" if suffixes_as_string else "Finding files"
    )
    logging.info(f"{file_message} corresponding to sample {sample} ...")

    if suffixes is None:
        # Return all files that contain `sample` in their name.
        matched_files = [f for f in files if sample in f.name]
        if not matched_files:
            raise ValueError(
                f"No files found containing sample '{sample}' with no suffix filtering."
            )
        return matched_files
    else:
        # We'll gather matches for each suffix separately, ensuring each suffix is found at least once.
        all_matched_files = []

        for suf in suffixes:
            # Collect files matching this particular suffix
            matched_for_suf = [
                f for f in files if sample in f.name and f.name.endswith(suf)
            ]
            if not matched_for_suf:
                # If *no* files for this suffix, raise an error immediately.
                raise ValueError(
                    f"No files found containing sample '{sample}' with suffix '{suf}'."
                )
            # Add these matches to the global list (we can deduplicate later if needed)
            all_matched_files.extend(matched_for_suf)

        return all_matched_files


# -----------------------------------------------
# Functions for processing the parameter file
# -----------------------------------------------
def load_metadata(parameter_json: Path) -> t.Dict[str, t.Any]:
    """Load JSON metadata from the given file."""
    logging.info(f"Loading metadata from {str(parameter_json)} ...")
    with open(parameter_json, "r") as json_file:
        return json.load(json_file)


def validate_metadata_keys(
    metadata: t.Dict[str, t.Any], expected_keys: t.Tuple[str, ...], parameter_json: Path
):
    """Validate that the expected keys are present in the metadata."""
    missing_keys = tuple(key for key in expected_keys if key not in metadata)
    if missing_keys:
        raise KeyError(
            f"The following keys are missing in {str(parameter_json)}: {', '.join(missing_keys)}. Please check input data."
        )


def convert_paths(
    metadata: t.Dict[str, t.Any], exclude_keys: t.Set[str]
) -> t.Dict[str, t.Any]:
    """Convert string file paths to Path objects, excluding specified keys."""
    metadata_with_paths = {}
    for key, value in metadata.items():
        if key not in exclude_keys:
            if isinstance(value, str):
                metadata_with_paths[key] = Path(value)
            elif isinstance(value, list) and all(isinstance(v, str) for v in value):
                metadata_with_paths[key] = [Path(v) for v in value]
            else:
                metadata_with_paths[
                    key
                ] = value  # Keep as-is for non-path-related values
        else:
            metadata_with_paths[key] = value  # Exclude from Path conversion
    return metadata_with_paths


def extract_metadata_files_from_parameter_json(
    parameter_json: Path,
) -> t.Dict[str, t.Any]:
    """Extract all metadata file paths from the parameter JSON and return as a dictionary."""
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
    exclude_keys = {"unplaced_contig_prefixes"}

    metadata = load_metadata(parameter_json)
    validate_metadata_keys(metadata, expected_keys, parameter_json)
    return convert_paths(metadata, exclude_keys)


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
    # Read in sample IDs from the exclude file
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


def expand_sex_abbreviation(sex: str) -> str:
    """
    Expands the sex abbreviations found in the metadata spreadsheet to their full word counterparts for better readability
    M -> male
    F -> female
    U -> unknown
    """
    if sex == "M":
        return "male"
    elif sex == "F":
        return "female"
    elif sex == "U":
        return "unknown"
    else:
        raise ValueError(f"Unrecognised value for sex: {sex}. Expected M, F or U.")


def get_sample_sex(sample_id: str, sample_metadata_xlsx: Path) -> str:
    """
    Extract the sex for a given sample ID from the metadata Excel file.

    Args:
        sample_id (str): The sample ID for which to retrieve the sex.
        sample_metadata_xlsx (Path): Path to the metadata Excel file.

    Returns:
        str: The sex of the given sample ID (expanded format: "Male", "Female", or "Unknown").

    Raises:
        ValueError: If the sample ID is not found or if the sex value is unexpected.
    """
    # Load the sample metadata spreadsheet
    sample_metadata_spreadsheet = pd.ExcelFile(sample_metadata_xlsx)

    # Iterate through each sheet (study) in the metadata spreadsheet
    for sheet_name in sample_metadata_spreadsheet.sheet_names:
        sheet_data = pd.read_excel(sample_metadata_spreadsheet, sheet_name=sheet_name)

        # Check if the sample ID exists in the current sheet
        sample_row = sheet_data[sheet_data["Sanger DNA ID"] == sample_id]
        if not sample_row.empty:
            # Extract the sex value
            sex = sample_row.iloc[0]["Sex"]
            if sex not in {"M", "F", "U"}:
                logger.error(
                    f"Unexpected sex value '{sex}' for sample '{sample_id}' in sheet '{sheet_name}'."
                )
                raise ValueError(
                    f"Unexpected sex value '{sex}' for sample '{sample_id}'. Expected one of 'M', 'F', or 'U'."
                )

            # Expand the sex abbreviation and return it
            return expand_sex_abbreviation(sex)

    # If the sample ID was not found in any sheet, raise an error
    logger.error(
        f"Sample ID '{sample_id}' is missing from the metadata Excel '{sample_metadata_xlsx}'."
    )
    raise ValueError(
        f"Sample ID '{sample_id}' is missing from the metadata Excel '{sample_metadata_xlsx}'."
    )


def determine_sample_sexes(
    file_list: t.List[Path], sample_metadata_xlsx: Path
) -> t.Dict[str, str]:
    """
    Determine the sex of each sample in the file list based on the metadata Excel file.
    """
    # Get the corresponding sample IDs for each file in the file list
    sample_ids = get_sample_ids_for_file_list(file_list)

    # Initialise a dictionary to store each sample's sex
    sample_sex_dict = dict()

    # Use the helper function to extract the sex for each sample ID
    for sample_id in sample_ids:
        try:
            sample_sex_dict[sample_id] = get_sample_sex(sample_id, sample_metadata_xlsx)
        except ValueError as e:
            logger.error(str(e))
            raise

    # Identify missing samples
    missing_samples = sample_ids - sample_sex_dict.keys()
    if missing_samples:
        missing_samples_str = ", ".join(sorted(missing_samples))
        logger.error(
            f"The following samples are missing from the metadata Excel '{sample_metadata_xlsx}': {missing_samples_str}."
        )
        raise ValueError(
            f"The following samples are missing from the metadata Excel '{sample_metadata_xlsx}': {missing_samples_str}."
        )

    return sample_sex_dict


def split_file_list_by_sample_sex(
    file_list: t.List[Path], sample_metadata_xlsx: Path
) -> t.DefaultDict[str, t.List[Path]]:
    """
    Splits a list of file paths by the sex of the samples they correspond to.
    This function takes a list of file paths and a metadata Excel file that contains
    information about the sex of each sample. It determines the sex of each sample
    and groups the file paths by sex.
    Args:
        file_list (t.List[Path]): A list of file paths to be split by sample sex.
        sample_metadata_xlsx (Path): Path to the Excel file containing sample metadata,
                                     including the sex of each sample.
    Returns:
        t.DefaultDict[str, t.List[Path]]: A dictionary where the keys are the sexes
                                          ('male', 'female', etc.) and the values are
                                          lists of file paths corresponding to samples
                                          of that sex.
    """
    # Determine the sexes of all samples in the file list
    sample_sex_dict = determine_sample_sexes(file_list, sample_metadata_xlsx)

    # Initialise a dictionary to store sex-separated files
    file_sex_dict = defaultdict(list)

    for file in file_list:
        sample_id = get_sample_id_from_file_path(file)
        sex = sample_sex_dict.get(sample_id)
        file_sex_dict[sex].append(file)

    return file_sex_dict


# -----------------------------------------------
# Functions for filtering by sample study
# -----------------------------------------------


def map_sample_ids_to_study_ids(
    sample_ids: set, sample_metadata_xlsx: Path
) -> t.Dict[str, str]:
    """
    Map sample IDs to their corresponding study IDs based on the metadata Excel file.

    Args:
        sample_ids (set): A set of sample IDs to map to study IDs.
        sample_metadata_xlsx (Path): Path to the metadata Excel file.

    Returns:
        Dict[str, str]: A dictionary mapping each sample ID to its corresponding study ID.
        Keys are study IDs and values are lists of sample IDs.

    Raises:
        ValueError: If a sample ID is not found in the metadata Excel file.
    """
    sample_study_dict = defaultdict(list)

    for sample_id in sample_ids:
        study_id = get_sample_study(sample_metadata_xlsx, sample_id)
        sample_study_dict[study_id].append(sample_id)

    return sample_study_dict


def get_sample_study(sample_metadata_xlsx: Path, sample_id: str) -> str:
    """Get the study ID for the given sample from the sample metadata Excel spreadsheet"""

    # Load the Excel file
    excel_data = pd.ExcelFile(sample_metadata_xlsx)

    # Iterate through each sheet in the file
    for sheet_name in excel_data.sheet_names:
        # Load the sheet into a DataFrame
        df = excel_data.parse(sheet_name)

        # Search for the sample ID
        matched_row = df[df["Sanger DNA ID"] == sample_id]

        # If a match is found, return the study ID
        if not matched_row.empty:
            return sheet_name

    # Raise an error if the sample ID is not found
    raise ValueError(f"Sample ID '{sample_id}' not found in any sheet.")


# -----------------------------------------------
# Functions for determining tumour/normal status of a given sample
# -----------------------------------------------
def get_tumour_normal_status(sample_metadata_xlsx: Path, sample_id: str):
    """
    This function takes an Excel file path and a sample ID, searches all sheets
    for the sample ID, and returns the tumour/normal status (T/N) of the sample.
    Raises a ValueError if the sample ID is not found, if required columns are missing,
    or if the T/N status is invalid.

    Parameters:
        sample_metadata_xlsx (Path): Path to the Excel file.
        sample_id (str): The sample ID to search for.

    Returns:
        str: Tumour/normal status ('T' or 'N') if found.

    Raises:
        ValueError: If the sample ID is not found in any sheet.
        ValueError: If required columns ('Sanger DNA ID', 'T/N') are missing in any sheet.
        ValueError: If the T/N status for the sample ID is invalid, i.e., does not start
                    with 'T' or 'N'.
    """
    # Load the Excel file
    excel_data = pd.ExcelFile(sample_metadata_xlsx)

    # Iterate through each sheet in the file
    for sheet_name in excel_data.sheet_names:
        # Load the sheet into a DataFrame
        df = excel_data.parse(sheet_name)

        # Check if the required columns exist
        if not {"Sanger DNA ID", "T/N"}.issubset(df.columns):
            raise ValueError(
                f"Required columns ('Sanger DNA ID', 'T/N') are missing in sheet '{sheet_name}'."
            )

        # Search for the sample ID
        matched_row = df[df["Sanger DNA ID"] == sample_id]

        # If a match is found, return the T/N status
        if not matched_row.empty:
            tn_status = matched_row["T/N"].iloc[0]
            if isinstance(tn_status, str) and tn_status.startswith(("T", "N")):
                return tn_status[
                    0
                ]  # Only return the first character (T or N) e.g. for samples with T1, N2 etc.
            else:
                raise ValueError(
                    f"Invalid T/N status '{tn_status}' for sample ID '{sample_id}'."
                )

    # Raise an error if the sample ID is not found
    raise ValueError(f"Sample ID '{sample_id}' not found in any sheet.")


# -----------------------------------------------
# Functions for categorising files based on their tumour/normal status
# -----------------------------------------------


def categorise_files_by_tumour_normal_status(
    files: t.List[Path], sample_metadata_xlsx: Path
) -> t.DefaultDict[str, t.List[Path]]:
    logging.info("Categorising files based on their tumour/normal status ...")

    # Initialise a defaultdict to store categorised files
    tn_status_file_dict = defaultdict(list)

    for file in files:
        sample_id = get_sample_id_from_file_path(file)
        tn_status = get_tumour_normal_status(sample_metadata_xlsx, sample_id)

        if tn_status == "T":
            tn_status_file_dict["T"].append(file)
        elif tn_status == "N":
            tn_status_file_dict["N"].append(file)
        else:
            logging.error(
                "Unable to categorise files based on their tumour/normal status"
            )
            raise ValueError(
                f"Got unexpected tumour normal status '{tn_status}' for sample ID '{sample_id}'. Expected 'T' or 'N'."
            )

    logging.debug(f"Categorised files: {tn_status_file_dict}")
    logging.info("Successfully categorised files by their tumour/normal status")

    return tn_status_file_dict
