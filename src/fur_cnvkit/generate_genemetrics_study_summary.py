#!/usr/bin/env python3
import argparse
from pathlib import Path
import typing as t

import pandas as pd

from fur_cnvkit import constants
from fur_cnvkit.utils.fur_utils import get_sample_id_from_file_path
from fur_cnvkit.utils.logging_utils import get_package_logger, setup_logging

COMMAND_NAME: str = constants.COMMAND_NAME__GENERATE_GENEMETRICS_STUDY_SUMMARY
logger = get_package_logger()


def get_argparser(
    subparser: t.Optional[argparse._SubParsersAction] = None,
) -> argparse.ArgumentParser:
    """
    Either returns a new ArgumentParser instance or a subparser for the
    genemetrics study summary command.
    """
    if subparser is None:
        parser = argparse.ArgumentParser(
            description=constants.DESCRIPTION__GENERATE_GENEMETRICS_STUDY_SUMMARY
        )
    else:
        parser = subparser.add_parser(
            COMMAND_NAME,
            description=constants.DESCRIPTION__GENERATE_GENEMETRICS_STUDY_SUMMARY,
            help=constants.SHORT_HELP__GENERATE_GENEMETRICS_STUDY_SUMMARY,
        )

    parser.add_argument(
        "-s",
        "--study-id",
        type=str,
        required=True,
        help="Study identifier used to prefix the output CSV.",
    )
    parser.add_argument(
        "-g",
        "--genemetrics-files",
        type=Path,
        nargs="+",
        default=None,
        help="Explicit genemetrics files to include (e.g. *.genemetrics.significant.out).",
    )
    parser.add_argument(
        "-d",
        "--genemetrics-dir",
        type=Path,
        default=None,
        help="Directory to search for genemetrics files (pattern controlled by --file-pattern).",
    )
    parser.add_argument(
        "--file-pattern",
        type=str,
        default="*.genemetrics.significant.out",
        help="Glob used when scanning --genemetrics-dir (default: %(default)s).",
    )
    parser.add_argument(
        "-b",
        "--baitset-genes",
        type=Path,
        required=True,
        help="Text file containing one baitset gene symbol per line.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for the study summary CSV.",
    )
    parser.add_argument(
        "--gain-threshold",
        type=float,
        default=None,
        help="Optional log2 threshold for gains to record in the summary.",
    )
    parser.add_argument(
        "--loss-threshold",
        type=float,
        default=None,
        help="Optional log2 threshold for losses to record in the summary.",
    )
    parser.add_argument(
        "--verbose",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        dest="verbose",
        help="log level",
    )
    return parser


def _resolve_explicit_genemetrics_files(
    explicit_files: t.Optional[t.Sequence[Path]],
) -> tuple[list[Path], list[Path]]:
    """
    Expand and validate explicit genemetrics file paths, returning valid and missing entries.
    """
    valid_files: list[Path] = []
    missing_files: list[Path] = []
    if not explicit_files:
        return valid_files, missing_files

    for file_path in explicit_files:
        expanded = file_path.expanduser()
        if expanded.is_file():
            valid_files.append(expanded)
        else:
            missing_files.append(expanded)
    return valid_files, missing_files


def _discover_genemetrics_files(
    directory: t.Optional[Path], pattern: str
) -> list[Path]:
    """
    Discover genemetrics files by globbing a directory with the provided pattern.
    """
    if directory is None:
        return []

    expanded_dir = directory.expanduser()
    if not expanded_dir.is_dir():
        raise ValueError(
            f"Genemetrics directory '{expanded_dir}' does not exist or is not a directory."
        )
    discovered = sorted(expanded_dir.glob(pattern))
    if not discovered:
        logger.warning(
            "Directory %s does not contain any files matching pattern '%s'.",
            expanded_dir,
            pattern,
        )
    return [path for path in discovered if path.is_file()]


def _deduplicate_paths(paths: t.Iterable[Path]) -> list[Path]:
    """
    Resolve duplicate file paths.
    """
    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        deduped.append(resolved)
    deduped.sort(key=get_sample_id_from_file_path)
    return deduped


def _filter_duplicate_samples(paths: t.Iterable[Path]) -> tuple[list[Path], list[str]]:
    """
    Ensure that only one genemetrics file per sample is returned, tracking duplicates.
    """
    unique_files: list[Path] = []
    duplicates: list[str] = []
    seen_samples: set[str] = set()

    for path in paths:
        sample_id = get_sample_id_from_file_path(path)
        if sample_id in seen_samples:
            duplicates.append(sample_id)
            continue
        seen_samples.add(sample_id)
        unique_files.append(path)

    return unique_files, duplicates


def _collect_genemetrics_files(
    explicit_files: t.Optional[t.Sequence[Path]],
    directory: t.Optional[Path],
    pattern: str,
) -> t.List[Path]:
    """
    Collect and de-duplicate genemetrics files supplied explicitly or via directory globbing.
    """
    explicit, missing = _resolve_explicit_genemetrics_files(explicit_files)
    discovered = _discover_genemetrics_files(directory, pattern)

    if missing:
        logger.warning(
            "Skipping %d missing genemetrics files: %s",
            len(missing),
            ", ".join(str(path) for path in missing),
        )

    all_candidates = explicit + discovered
    if not all_candidates:
        raise ValueError(
            "No genemetrics files found. Provide --genemetrics-files or --genemetrics-dir with matching files."
        )

    deduped = _deduplicate_paths(all_candidates)
    unique_files, duplicates = _filter_duplicate_samples(deduped)

    if duplicates:
        logger.warning(
            "Multiple genemetrics files detected for samples: %s. Using the first occurrence of each.",
            ", ".join(sorted(set(duplicates))),
        )

    logger.info("Using %d genemetrics files for summary generation.", len(unique_files))
    return unique_files


def _read_baitset_genes(baitset_genes_file: Path) -> list[str]:
    """
    Return a list of baitset genes from the provided file, ensuring it is not empty.
    """
    with baitset_genes_file.open("r") as handle:
        genes = [line.strip() for line in handle if line.strip()]
    if not genes:
        raise ValueError(
            f"No gene symbols were found in baitset file '{baitset_genes_file}'."
        )
    return genes


def _load_sample_series(
    genemetrics_file: Path, baitset_genes: t.Sequence[str], skip_labels: set[str]
) -> pd.Series:
    """
    Load gene-level log2 values for a single sample into a Series indexed by baitset genes.
    """
    sample_id = get_sample_id_from_file_path(genemetrics_file)
    logger.info("Processing genemetrics file for sample %s ...", sample_id)
    df_temp = pd.read_csv(genemetrics_file, sep="\t", usecols=[0, 4])
    logger.debug("Initial data for sample %s:\n%s", sample_id, df_temp.head())

    df_temp = df_temp[~df_temp["gene"].isin(skip_labels)]
    df_temp["log2"] = pd.to_numeric(df_temp["log2"], errors="coerce")
    df_temp = df_temp.groupby("gene", as_index=False).mean()
    df_temp.set_index("gene", inplace=True)
    df_temp = df_temp.reindex(baitset_genes)
    return df_temp["log2"].rename(sample_id)


def _apply_threshold_columns(
    cohort_df: pd.DataFrame,
    gain_threshold: t.Optional[float],
    loss_threshold: t.Optional[float],
) -> pd.DataFrame:
    """
    Add gain/loss threshold columns when thresholds are provided.
    """
    if gain_threshold is not None:
        cohort_df["gain_threshold"] = gain_threshold
    if loss_threshold is not None:
        cohort_df["loss_threshold"] = loss_threshold
    return cohort_df


def generate_genemetrics_study_summary_csv(
    study_id: str,
    genemetrics_files: t.Sequence[Path],
    baitset_genes_file: Path,
    outdir: Path,
    gain_threshold: t.Optional[float] = None,
    loss_threshold: t.Optional[float] = None,
) -> pd.DataFrame:
    """
    Generate a study-level CSV summarizing genemetrics data.

    The resulting CSV contains one row per sample, gene-level log2 values as columns,
    and optional `gain_threshold`/`loss_threshold` columns when thresholds are supplied.
    """
    logger.info("Generating study summary CSV for study %s ...", study_id)
    outdir.mkdir(parents=True, exist_ok=True)
    output_csv_path = outdir / f"{study_id}.genemetrics_study_summary.csv"
    logger.debug("Output CSV path: %s", output_csv_path)

    baitset_genes_list = _read_baitset_genes(baitset_genes_file)
    logger.debug("Loaded %d baitset genes.", len(baitset_genes_list))

    sample_series_list: list[pd.Series] = []
    skip_labels = {"none", "unk", "incmpl", "cmpl"}

    for genemetrics_file in genemetrics_files:
        sample_series_list.append(
            _load_sample_series(genemetrics_file, baitset_genes_list, skip_labels)
        )

    if not sample_series_list:
        raise ValueError(
            "No usable genemetrics data found in the supplied files after filtering."
        )

    cohort_df = pd.concat(sample_series_list, axis=1).T
    cohort_df = _apply_threshold_columns(
        cohort_df, gain_threshold=gain_threshold, loss_threshold=loss_threshold
    )
    cohort_df.to_csv(output_csv_path)
    logger.info("Study summary CSV created at %s", output_csv_path)
    return cohort_df


def main(args: t.Optional[argparse.Namespace] = None) -> t.Optional[int]:
    """
    Entry point for CLI execution.
    """
    if args is None:
        args = get_argparser().parse_args()
    logger.debug("Parsed arguments: %s", args)

    genemetrics_files = _collect_genemetrics_files(
        args.genemetrics_files,
        args.genemetrics_dir,
        args.file_pattern,
    )
    generate_genemetrics_study_summary_csv(
        study_id=args.study_id,
        genemetrics_files=genemetrics_files,
        baitset_genes_file=args.baitset_genes,
        outdir=args.outdir,
        gain_threshold=args.gain_threshold,
        loss_threshold=args.loss_threshold,
    )
    return None


if __name__ == "__main__":
    parsed_args = get_argparser().parse_args()
    setup_logging(level=getattr(parsed_args, "verbose", "INFO"))
    main(parsed_args)
