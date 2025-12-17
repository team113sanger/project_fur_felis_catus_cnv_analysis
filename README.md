# fur_cnvkit

This repository contains code that was used for performing a copy number alteration (CNA) analysis for the FUR Felis catus project in the Adams group at the Wellcome Trust Sanger institute. It builds on top of [CNVkit](http://dx.doi.org/10.1371/journal.pcbi.1004873) and can be adapted for other targeted sequencing projects that require a similar analysis by altering the parameters in `src/fur_cnvkit/parameters/parameter_file_template.json`

|                         Main                         |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][main-pipe-badge]][main-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[main-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_cnvkit/badges/main/pipeline.svg
[main-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_cnvkit/-/commits/main
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_cnvkit/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_cnvkit/-/commits/develop


## Usage

Within the Docker image (see below) you can run the following command:

```
usage: fur_cnvkit [-h] [--version] COMMAND ...

FUR CNVkit - a collection commands for CN analysis as part of FUR project that
rely on CNVkit.

positional arguments:
  COMMAND
    mad                 Calculate the weighted Median Absolute Deviation (MAD) of segments in samples.
    cnvkit_static_files
                        Generate CNVKit static files.
    copy_number_reference
                        Generate a copy number reference file.
    oncoprint           Generate an oncoprint figure.
    cnvkit_cn_calling_pipeline
                        Run the end-to-end copy number calling workflow for a cohort of samples

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

Please see tests for examples

### CPU planning for pipeline runs

Several commands expose knobs that let you balance throughput against overall CPU usage:

- `generate_copy_number_reference --threads N` limits the number of worker processes used when building coverage files and running normal-vs-normal comparisons. Setting it to the total number of CPUs on the host (e.g. `--threads 32`) gives CNVkit free rein; lower values reduce concurrency.
- `run_cnvkit_copy_number_calling_pipeline --study-workers S --batch-processes B` lets you bound the cross-study fan-out (`S`) and how many processes each `cnvkit.py batch` invocation uses internally (`B`). The upper bound on simultaneous CNVkit workers is roughly `S × B`.

Practical guidance:

1. Decide how many studies can run side-by-side without thrashing the filesystem or saturating memory, and pass that to `--study-workers`. Leaving it unset keeps the legacy sequential behaviour.
2. For each study, pick a `--batch-processes` value that keeps `S × B` below the total CPU count (e.g. on a 48‑core host, `--study-workers 3 --batch-processes 12` consumes at most 36 cores, leaving headroom for the OS and I/O threads).
3. If you want to devote the entire machine to a single study, set `--study-workers 1 --batch-processes <CPU count>` (or omit `--batch-processes` to let CNVkit auto-detect and use all CPUs).

Adjust these numbers downward if you observe I/O bottlenecks or if other services share the machine.

## Docker Image

This project hosts Docker images on Quay.io. Please see [https://quay.io/repository/team113sanger/fur_cnvkit](https://quay.io/repository/team113sanger/fur_cnvkit?tab=tags).

## Summary

Given
- A cohort of tumor and normal sample `.bam` files and their indexes
- A bedfile describing the targeted regions in a pulldown study
- A reference genome fasta file
- A flat file containing gene coordinates (e.g. Ensembl GTF or GFF3)
- A sample metadata table (Excel/TSV/CSV) containing sample IDs (`Sanger DNA ID` by default), tumour/normal status (`T/N`) and Sex (`Sex`). Column names can be customised via the CNVKit static file CLI or by setting the `metadata_columns` entry in the generated parameter file if your manifests use different headers (for example `Sanger_DNA_ID` or `Phenotype`).


This tool:
- Splits the cohort into male and female sub-cohorts
- Performs `cnvkit assess` and `cnvkit autobin` to determine mappable regions of the genome and split the genome into target and anti-target regions (generating static files).
- Determines which normal samples in each sub-cohort have copy number profiles of sufficient quality to include in a pooled copy number reference by normal vs. normal copy number calling ([Chandramohan et al, 2022](https://pubmed.ncbi.nlm.nih.gov/35487348/)). Samples are evaluated based on the difference between the **median gene-level log2 ratio** in a sample's copy number profile vs the pooled reference of all other samples and the **median of the median gene-level log2 ratios** for all samples in the sub-cohort.
- Creates the pooled copy-number reference file for a sub-cohort from the highest quality normals, discarding the bottom 20% of normal samples.
- Calls copy number alterations in tumour samples against the pooled reference using the `cnvikit batch` sub-command
- Calculates the median absolute deviation (MAD) of the log₂ copy number ratios from segments in each sample, weighted the by average segment size to identify tumours with low-quality, hyper-segmented copy number profiles.
- Removes tumours with low-quality and hypersegmented copy number profiles from each subcohort and optionally generates an oncoprint figure for the remaining samples


## Directory structure

```
.
├── .devcontainer
│   └── devcontainer.json    # Used by VSCode to create a development container
├── cicd_scripts/
│   ├── build_package_and_publish.sh       # Publishes the package to the Sanger Gitlab PyPi registry
│   └── format_test.sh       # CI/CD script to run formatting and linting
├── .dockerignore
├── .gitignore
├── .gitlab-ci.yml           # Gitlab CI/CD config file
├── .pre-commit-config.yaml  # Pre-commit config file
├── Dockerfile
├── README.md
├── docker-compose.yml
├── pyproject.toml           # Python package config file (dependencies add here)
├── poetry.lock              # Poetry dependency lock file (generated)
├── poetry.toml              # Poetry config file
├── requirements.txt         # Alternative package dependency file (generated)
├── src
│   └── fur_cnvkit  # The python package of this repo (src code)
└── tests
    ├── integration          # Integration tests go here
    └── unit                 # Unit tests go here
```

## CI/CD Pipeline Overview

The project's CI/CD pipeline is configured in `.gitlab-ci.yml` and comprises three stages: `build`, `test`, and `publish`.

### Local/Remote Machine Setup

#### Prerequisites:
 - Python 3.10 (or later) installed
 - Poetry installed
    - Follow the [official instructions here](https://python-poetry.org/docs/#installation)
- Git + hubflow installed

## Running Tests

This project uses pytest as the main testing framework and all tests are located in the `tests` directory. For formatting and linting, we use `pre-commit` hooks which are installed in the setup steps above.

### Unit & Integration Tests

Before running the tests, ensure that your virtual environment is activated and all the required packages are installed. If you followed the above setup steps, your environment should be ready.

To run the tests, navigate to the project's root directory and execute the following command:

```bash
# To run all tests
pytest
# If you want to see a more verbose output
pytest -sv
# If you want to have the tests stop after the second failure
pytest --maxfail 2
# If you want to run only the failing tests from the last run
pytest --lf
```

### Pre-commit (Linting & Formatting)

This project uses `pre-commit` hooks to ensure that all code is formatted and linted before committing. That means whenever you commit the code will be checked and then check again on the remote GitLab server. If the code is not formatted or linted, the commit will fail.

To run the hooks manually, execute the following command:

```bash
pre-commit run --all-files
# Or use the CICD compatible script
./cicd_scripts/format_test.sh
```

## Rolling out a Tagged Release
This repo the GitFlow branching model and uses [hubflow](https://datasift.github.io/gitflow/TheHubFlowTools.html) as a tool to enable this from the CLI.

```bash
# Switch to the develop branch
git checkout develop

# Start a new release branch e.g. 0.1.0 not v0.1.0
git hf release start <project_version>
```

Now, do the following things:
* `CHANGELOG.md`: Under the heading of the newest release version, describe what was changed, fixed, added.
* `pyproject.toml`: Increment the project version to the current release version
* Commit these changes

Finally

```bash
git hf release finish <project_version>
```
