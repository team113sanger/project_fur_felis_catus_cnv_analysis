# This is where you may define fixtures that are shared across multiple test files.
#
# See this FAQ: https://docs.pytest.org/en/8.2.x/how-to/index.html
# Or read this section of the docs: https://docs.pytest.org/en/8.2.x/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import logging
import os
import pathlib

import pytest


# CONSTANTS
ENV_VAR_BAM_DIR = "TEST_BAM_DIR"
ENV_VAR_BAITSET_DIR = "TEST_BAITSET_DIR"
ENV_VAR_GENOME_DIR = "TEST_GENOME_DIR"
ENV_VAR_ANNOTATION_DIR = "TEST_ANNOTATION_DIR"
ENV_VAR_CNVKIT_METADATA_DIR = "TEST_CNVKIT_METADATA_DIR"


# FIXTURES
@pytest.fixture
def mock_file_path(tmp_path):
    """Fixture to create a temporary file for testing."""
    return tmp_path / "test_file"


@pytest.fixture
def bams() -> list[pathlib.Path]:
    raw_bam_dir = os.environ.get(ENV_VAR_BAM_DIR, "")
    bam_dir = pathlib.Path(raw_bam_dir)

    return list(bam_dir.glob("**/*.bam"))


@pytest.fixture
def feline_baitset() -> pathlib.Path:
    raw_baitset_dir = os.environ.get(ENV_VAR_BAITSET_DIR, "")
    expected_baitset = "S3250994_Feline_HSA_Jan2020_146_canonical_pad100.merged.bed"
    baitset = pathlib.Path(raw_baitset_dir) / expected_baitset
    if not baitset.exists():
        raise FileNotFoundError(f"Could not find baitset at {baitset}")
    return baitset


@pytest.fixture
def feline_reference_fasta() -> pathlib.Path:
    raw_genome_dir = os.environ.get(ENV_VAR_GENOME_DIR, "")
    expected_reference_fasta = "Felis_catus.Felis_catus_9.0.dna.toplevel.fa"
    reference_fasta = pathlib.Path(raw_genome_dir).resolve() / expected_reference_fasta
    if not reference_fasta.exists():
        raise FileNotFoundError(
            f"Could not find reference FASTA at {str(reference_fasta)}"
        )
    return reference_fasta


@pytest.fixture
def feline_refflat_file() -> pathlib.Path:
    raw_annotation_dir = os.environ.get(ENV_VAR_ANNOTATION_DIR, "")
    if not raw_annotation_dir:
        raise ValueError(
            f"Environment variable {ENV_VAR_ANNOTATION_DIR} is not set or empty."
        )
    refflat_dir = pathlib.Path(raw_annotation_dir).resolve() / "refFlat_files"
    expected_refflat_file = "Felis_catus.Felis_catus_9.0.104.refFlat.txt"
    refflat_file = refflat_dir / expected_refflat_file
    if not refflat_file.exists():
        raise FileNotFoundError(f"Could not find refFlat file at {str(refflat_file)}")
    return refflat_file


@pytest.fixture
def feline_cnvkit_access_bed() -> pathlib.Path:
    raw_cnvkit_metadata_dir = os.environ.get(ENV_VAR_CNVKIT_METADATA_DIR, "")
    cnvkit_access_dir = pathlib.Path(raw_cnvkit_metadata_dir).resolve() / "access"
    expected_cnvkit_access_bed = "access-Felis_catus.Felis_catus_9.0.dna.toplevel.bed"
    cnvkit_access_bed = cnvkit_access_dir / expected_cnvkit_access_bed
    if not cnvkit_access_bed.exists():
        raise FileNotFoundError(
            f"Could not find CNVKit access BED at {str(cnvkit_access_bed)}"
        )
    return cnvkit_access_bed
