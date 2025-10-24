from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
import typing as t

import pandas as pd

from fur_cnvkit.utils.logging_utils import get_package_logger

logger = get_package_logger()


# -----------------------------------------------
# Metadata configuration handling
# -----------------------------------------------
@dataclass(frozen=True)
class MetadataColumns:
    """Container for the column names used in sample metadata files."""

    sample_id: str = "Sanger DNA ID"
    tumour_normal: str = "T/N"
    sex: str = "Sex"
    study: t.Optional[str] = None


DEFAULT_METADATA_COLUMNS = MetadataColumns()
_metadata_columns: MetadataColumns = DEFAULT_METADATA_COLUMNS
_metadata_cache: t.Dict[t.Tuple[Path, MetadataColumns], t.Dict[str, pd.DataFrame]] = {}


def _normalise_metadata_config(
    config: t.Optional[t.Union[MetadataColumns, t.Dict[str, t.Optional[str]]]]
) -> MetadataColumns:
    """Normalise user supplied metadata column configuration."""
    if config is None:
        return DEFAULT_METADATA_COLUMNS
    if isinstance(config, MetadataColumns):
        return config

    if not isinstance(config, dict):
        raise TypeError(
            "Metadata column configuration must be a MetadataColumns instance, a dict, or None."
        )

    allowed_keys = {"sample_id", "tumour_normal", "sex", "study"}
    invalid_keys = set(config.keys()) - allowed_keys
    if invalid_keys:
        invalid = ", ".join(sorted(invalid_keys))
        raise KeyError(
            f"Unknown metadata column keys provided: {invalid}. "
            f"Supported keys are: {', '.join(sorted(allowed_keys))}."
        )

    return MetadataColumns(
        sample_id=config.get("sample_id", DEFAULT_METADATA_COLUMNS.sample_id),
        tumour_normal=config.get(
            "tumour_normal", DEFAULT_METADATA_COLUMNS.tumour_normal
        ),
        sex=config.get("sex", DEFAULT_METADATA_COLUMNS.sex),
        study=config.get("study", DEFAULT_METADATA_COLUMNS.study),
    )


def set_metadata_columns(
    config: t.Optional[t.Union[MetadataColumns, t.Dict[str, t.Optional[str]]]]
) -> MetadataColumns:
    """
    Configure the column names used when parsing sample metadata.

    Args:
        config: Either a MetadataColumns instance, a dictionary containing any of
            the keys 'sample_id', 'tumour_normal', 'sex', 'study', or None to reset
            to the defaults.

    Returns:
        The effective MetadataColumns configuration.
    """
    global _metadata_columns
    _metadata_columns = _normalise_metadata_config(config)
    _metadata_cache.clear()
    return _metadata_columns


def get_metadata_columns() -> MetadataColumns:
    """Return the current metadata column configuration."""
    return _metadata_columns


# -----------------------------------------------
# Metadata loading helpers
# -----------------------------------------------


def _ensure_required_columns(
    df: pd.DataFrame,
    required_columns: t.Iterable[str],
    metadata_path: Path,
    context: str,
) -> None:
    """Verify that the required columns exist in the provided DataFrame."""
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        missing_cols = ", ".join(missing)
        raise ValueError(
            f"Required column(s) {missing_cols} are missing from metadata source "
            f"'{metadata_path}' (context: {context})."
        )


def _read_tabular_file(metadata_path: Path) -> pd.DataFrame:
    """Load a tabular (TSV/CSV) metadata file."""
    suffix = metadata_path.suffix.lower()
    if suffix in (".tsv", ".txt"):
        return pd.read_csv(metadata_path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(metadata_path)

    raise ValueError(
        f"Unsupported metadata file extension '{suffix}' for '{metadata_path}'. "
        "Supported formats: .xlsx, .xls, .tsv, .txt, .csv."
    )


def _group_dataframe_by_study(
    df: pd.DataFrame, metadata_path: Path
) -> t.Dict[str, pd.DataFrame]:
    """Split a DataFrame by the configured study column."""
    columns = get_metadata_columns()
    if columns.study is None:
        return {"default": df.copy()}

    study_column = columns.study
    _ensure_required_columns(df, [study_column], metadata_path, "group by study")
    grouped: t.Dict[str, pd.DataFrame] = {}
    for study_id, group_df in df.groupby(study_column):
        grouped[str(study_id)] = group_df.copy()
    return grouped


def _load_metadata_tables(metadata_path: Path) -> t.Dict[str, pd.DataFrame]:
    """Load a metadata file into a dictionary of DataFrames keyed by study."""
    columns = get_metadata_columns()
    resolved_path = metadata_path.resolve()

    if not resolved_path.exists():
        raise FileNotFoundError(
            f"Metadata file '{resolved_path}' does not exist. Please check the path."
        )

    suffix = resolved_path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        spreadsheet = pd.ExcelFile(resolved_path)
        sheet_tables = {
            sheet_name: spreadsheet.parse(sheet_name).copy()
            for sheet_name in spreadsheet.sheet_names
        }
        if columns.study is not None:
            combined_df = pd.concat(sheet_tables.values(), ignore_index=True)
            return _group_dataframe_by_study(combined_df, resolved_path)
        return sheet_tables

    tabular_df = _read_tabular_file(resolved_path)
    return _group_dataframe_by_study(tabular_df, resolved_path)


def _get_metadata_tables(metadata_path: Path) -> t.Dict[str, pd.DataFrame]:
    """
    Return cached metadata tables keyed by study.

    Returns copies of the cached DataFrames to avoid accidental mutation.
    """
    cache_key = (metadata_path.resolve(), get_metadata_columns())
    if cache_key not in _metadata_cache:
        _metadata_cache[cache_key] = _load_metadata_tables(metadata_path)

    return {study: df.copy() for study, df in _metadata_cache[cache_key].items()}


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
    logger.info(f"{file_message} corresponding to sample {sample} ...")

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
    logger.info(f"Loading metadata from {str(parameter_json)} ...")
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
    if "sample_metadata_xlsx" not in metadata:
        if "sample_metadata_file" in metadata:
            metadata["sample_metadata_xlsx"] = metadata["sample_metadata_file"]
        else:
            raise KeyError(
                f"Parameter file '{parameter_json}' must contain either "
                "'sample_metadata_xlsx' or 'sample_metadata_file'."
            )
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
        logger.error(
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
    logger.info(f"Extracting samples to exclude from {str(exclude_file)} ...")
    exclude_samples = extract_sample_ids_from_exclude_file(exclude_file)

    logger.debug(f"Exclude samples: {exclude_samples}")

    # Filter out any files that belong to unwanted samples
    logger.info("Filtering out samples found in the exclude file ...")
    filtered_files = filter_files_by_exclude_samples(files, exclude_samples)

    logger.info(
        f"Successfully filtered out {len(files) - len(filtered_files)} files from unwanted samples, leaving a total of {len(filtered_files)} files."
    )
    logger.debug(f"Filtered files: {filtered_files}")

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


def get_sample_sex(sample_id: str, sample_metadata_file: Path) -> str:
    """
    Extract the sex for a given sample ID from the metadata file.

    Args:
        sample_id (str): The sample ID for which to retrieve the sex.
        sample_metadata_file (Path): Path to the metadata file.

    Returns:
        str: The sex of the given sample ID (expanded format: "male", "female", or "unknown").

    Raises:
        ValueError: If the sample ID is not found or if the sex value is unexpected.
    """
    columns = get_metadata_columns()
    tables = _get_metadata_tables(sample_metadata_file)

    for context, table in tables.items():
        _ensure_required_columns(
            table,
            (columns.sample_id, columns.sex),
            sample_metadata_file,
            context,
        )
        sample_row = table[table[columns.sample_id] == sample_id]
        if not sample_row.empty:
            sex_value = sample_row.iloc[0][columns.sex]
            if sex_value not in {"M", "F", "U"}:
                logger.error(
                    f"Unexpected sex value '{sex_value}' for sample '{sample_id}' in '{sample_metadata_file}' (context: {context})."
                )
                raise ValueError(
                    f"Unexpected sex value '{sex_value}' for sample '{sample_id}'. Expected one of 'M', 'F', or 'U'."
                )
            return expand_sex_abbreviation(sex_value)

    logger.error(
        f"Sample ID '{sample_id}' is missing from the metadata file '{sample_metadata_file}'."
    )
    raise ValueError(
        f"Sample ID '{sample_id}' is missing from the metadata file '{sample_metadata_file}'."
    )


def determine_sample_sexes(
    file_list: t.List[Path], sample_metadata_file: Path
) -> t.Dict[str, str]:
    """
    Determine the sex of each sample in the file list based on the metadata file.
    """
    # Get the corresponding sample IDs for each file in the file list
    sample_ids = get_sample_ids_for_file_list(file_list)

    # Initialise a dictionary to store each sample's sex
    sample_sex_dict = dict()

    # Use the helper function to extract the sex for each sample ID
    for sample_id in sample_ids:
        try:
            sample_sex_dict[sample_id] = get_sample_sex(sample_id, sample_metadata_file)
        except ValueError as e:
            logger.error(str(e))
            raise

    # Identify missing samples
    missing_samples = sample_ids - sample_sex_dict.keys()
    if missing_samples:
        missing_samples_str = ", ".join(sorted(missing_samples))
        logger.error(
            f"The following samples are missing from the metadata file '{sample_metadata_file}': {missing_samples_str}."
        )
        raise ValueError(
            f"The following samples are missing from the metadata file '{sample_metadata_file}': {missing_samples_str}."
        )

    return sample_sex_dict


def split_file_list_by_sample_sex(
    file_list: t.List[Path], sample_metadata_file: Path
) -> t.DefaultDict[str, t.List[Path]]:
    """
    Splits a list of file paths by the sex of the samples they correspond to.
    This function takes a list of file paths and a metadata file that contains
    information about the sex of each sample. It determines the sex of each sample
    and groups the file paths by sex.
    Args:
        file_list (t.List[Path]): A list of file paths to be split by sample sex.
        sample_metadata_file (Path): Path to the file containing sample metadata,
                                     including the sex of each sample.
    Returns:
        t.DefaultDict[str, t.List[Path]]: A dictionary where the keys are the sexes
                                          ('male', 'female', etc.) and the values are
                                          lists of file paths corresponding to samples
                                          of that sex.
    """
    # Determine the sexes of all samples in the file list
    sample_sex_dict = determine_sample_sexes(file_list, sample_metadata_file)

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
    sample_ids: set, sample_metadata_file: Path
) -> t.Dict[str, str]:
    """
    Map sample IDs to their corresponding study IDs based on the metadata file.

    Args:
        sample_ids (set): A set of sample IDs to map to study IDs.
        sample_metadata_file (Path): Path to the metadata file.

    Returns:
        Dict[str, str]: A dictionary mapping each sample ID to its corresponding study ID.
        Keys are study IDs and values are lists of sample IDs.

    Raises:
        ValueError: If a sample ID is not found in the metadata Excel file.
    """
    sample_study_dict = defaultdict(list)

    for sample_id in sample_ids:
        study_id = get_sample_study(sample_metadata_file, sample_id)
        sample_study_dict[study_id].append(sample_id)

    return sample_study_dict


def get_sample_study(sample_metadata_file: Path, sample_id: str) -> str:
    """Get the study ID for the given sample from the metadata file."""
    columns = get_metadata_columns()
    tables = _get_metadata_tables(sample_metadata_file)

    for context, table in tables.items():
        _ensure_required_columns(
            table,
            (columns.sample_id,),
            sample_metadata_file,
            context,
        )
        matched_row = table[table[columns.sample_id] == sample_id]

        if not matched_row.empty:
            if columns.study and columns.study in matched_row.columns:
                study_value = matched_row.iloc[0][columns.study]
                return str(study_value)
            return context

    raise ValueError(f"Sample ID '{sample_id}' not found in metadata.")


# -----------------------------------------------
# Functions for determining tumour/normal status of a given sample
# -----------------------------------------------
def get_tumour_normal_status(sample_metadata_file: Path, sample_id: str):
    """
    This function takes a metadata file path and a sample ID, searches all studies
    for the sample ID, and returns the tumour/normal status of the sample. Raises a
    ValueError if the sample ID is not found, if required columns are missing,
    or if the tumour/normal status is invalid.

    Parameters:
        sample_metadata_file (Path): Path to the metadata file.
        sample_id (str): The sample ID to search for.

    Returns:
        str: Tumour/normal status ('T' or 'N') if found.

    Raises:
        ValueError: If the sample ID is not found in any context.
        ValueError: If required columns are missing in any context.
        ValueError: If the T/N status for the sample ID is invalid, i.e., does not start
                    with 'T' or 'N'.
    """
    columns = get_metadata_columns()
    tables = _get_metadata_tables(sample_metadata_file)

    for context, table in tables.items():
        _ensure_required_columns(
            table,
            (columns.sample_id, columns.tumour_normal),
            sample_metadata_file,
            context,
        )

        matched_row = table[table[columns.sample_id] == sample_id]

        if not matched_row.empty:
            tn_status = matched_row.iloc[0][columns.tumour_normal]
            if isinstance(tn_status, str) and tn_status.startswith(("T", "N")):
                return tn_status[0]
            raise ValueError(
                f"Invalid {columns.tumour_normal} status '{tn_status}' for sample ID '{sample_id}'."
            )

    raise ValueError(f"Sample ID '{sample_id}' not found in metadata.")


# -----------------------------------------------
# Functions for categorising files based on their tumour/normal status
# -----------------------------------------------


def categorise_files_by_tumour_normal_status(
    files: t.List[Path], sample_metadata_file: Path
) -> t.DefaultDict[str, t.List[Path]]:
    logger.info("Categorising files based on their tumour/normal status ...")

    # Initialise a defaultdict to store categorised files
    tn_status_file_dict = defaultdict(list)

    for file in files:
        sample_id = get_sample_id_from_file_path(file)
        tn_status = get_tumour_normal_status(sample_metadata_file, sample_id)

        if tn_status == "T":
            tn_status_file_dict["T"].append(file)
        elif tn_status == "N":
            tn_status_file_dict["N"].append(file)
        else:
            logger.error(
                "Unable to categorise files based on their tumour/normal status"
            )
            raise ValueError(
                f"Got unexpected tumour normal status '{tn_status}' for sample ID '{sample_id}'. Expected 'T' or 'N'."
            )

    logger.debug(f"Categorised files: {tn_status_file_dict}")
    logger.info("Successfully categorised files by their tumour/normal status")

    return tn_status_file_dict
