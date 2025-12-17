# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2025-12-17
### Added
- `run_cnvkit_copy_number_calling_pipeline` can now process studies in parallel via `--study-workers` and limit CNVkit's internal multiprocessing per study with `--batch-processes`.

## [0.5.0] - 2025-12-10
### Added
- `generate_cnvkit_static_files` can now run in `--parameter-file-only` mode with new `--existing-*` arguments, allowing reuse of pre-generated BED and baitset files when only refreshing the `parameters.json`.
- `run_cnvkit_copy_number_calling_pipeline` gained `--weight_filter_threshold`, which filters low-weight segments from median-centred `.cns` files and emits tagged scatter/diagram plots built from the filtered calls.
- Introduced the reusable `filter_cns_by_weight` helper (with unit tests) plus support for tagged plot filenames in `run_cnvkit_diagram` / `run_cnvkit_scatter` to underpin the new pipeline option.

## [0.4.0] - 2025-11-17
### Added
- Added compatability with DERMATLAS-style metadata manifest files in `generate_cnvkit_static_files` instead of FUR-based Excel spreadsheets
- The `run_cnvkit_copy_number_calling_pipeline` step will now generate an unfiltered genemetrics summary file that no longer exludes outlier MAD-flagged samples
- Introduced a standalone `generate_genemetrics_study_summary` CLI that can recreate the per-study cohort summaries outside the full CNVKit pipeline.
- `--threads` option added to `generate_copy_number_reference` to allow parallelisation and speed up pipeline

## [0.3.0] - 2025-10-08
### Added
- CLI argument `generate_cnvkit_static_files -x` for passing a bedfile through to `cnvkit.py access -x` and refining the set of antitarget regions generated

## [0.2.0] - 2025-04-29
### Added
- Added a unified CLI with a single entrypoint and improved command structure across all scripts
- Added oncoprint visualization with sample filtering, clustering capabilities, and gene label formatting
- Added copy number reference generation pipeline with normal vs. normal comparisons
- Added gene-to-chrom mapping file and sample filtering by MAD (Median Absolute Deviation)
- Added functionality to determine and reclassify sample sex from "unknown" to "male" or "female"
- Added generation of cohort-level genemetrics summary file
- Added plot sorting options by number of alterations or genomic position
- Added R dependencies including ComplexHeatmap and integration tests

### Changed
- Improved logging system with verbose options and package-wide logger configuration
- Enhanced parallel processing with `max_cpus` argument for CNVKit subprocessing operations
- Swapped mode centering for median centering in call.cns files
- Improved sample metadata handling and parameter file generation
- Enhanced R environment handling with better package installation reliability
- Updated from Python 3.9 to Python 3.10

### Fixed
- Fixed parameter file naming conventions and argument handling
- Fixed issues with subcommands using sub-argparsers
- Fixed failing tests related to expected prefix JSON naming convention

## [0.1.0] - 2024-11-01
### Added
- Initial project setup using the template from [example-python-cicd](https://gitlab.internal.sanger.ac.uk/team113sanger/common/example-python-cicd).
- Template version is 0.1.4 -- consult its `CHANGELOG.md` for more details.
