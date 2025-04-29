# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
