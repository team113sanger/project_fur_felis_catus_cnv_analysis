import shutil
import pathlib
import os
import subprocess

import pytest

# CONSTANTS
ENV_VAR_BAM_DIR = "TEST_BAM_DIR"
ENV_VAR_BAITSET_DIR = "TEST_BAITSET_DIR"

# HELPERS


def can_find_bams() -> bool:
    raw_bam_dir = os.environ.get(ENV_VAR_BAM_DIR, "")
    bam_dir = pathlib.Path(raw_bam_dir)

    return bool(raw_bam_dir) and bam_dir.exists() and any(bam_dir.glob("**/*.bam"))


def can_find_baitset() -> bool:
    raw_baitset = os.environ.get(ENV_VAR_BAITSET_DIR, "")
    baitset = pathlib.Path(raw_baitset)

    return bool(raw_baitset) and baitset.exists() and any(baitset.glob("**/*.bed"))


def should_skip_tests() -> bool:
    return not all([can_find_bams(), can_find_baitset()])


# FIXTURES


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


# TESTS


def test_cnvkit_is_installed():
    program = "cnvkit.py"
    assert bool(shutil.which(program))


def test_cnvkit_version_is_expected():
    # Given
    expected_version = "0.9.10"
    cmd = "cnvkit.py version"

    # When
    version = subprocess.check_output(cmd, shell=True).decode().strip()

    # Then
    assert version == expected_version


@pytest.mark.skipif(should_skip_tests(), reason="No test data found")
def test_cnvkit_runs__numpy_regression():
    # Given
    pass  ## WIP TEST

    # When

    # Then


@pytest.mark.skipif(should_skip_tests(), reason="No test data found")
def test_cnvkit_runs__target(feline_baitset: pathlib.Path, tmp_path: pathlib.Path):
    # Given
    output_file = tmp_path / "feline_targets.bed"
    cmd = f"cnvkit.py target {feline_baitset} -o {str(output_file)}"

    expected_first_line = "X\t912693\t914028\t-"
    expected_stderr_substrings = [
        "Detected file format: bed",
        f"{str(output_file)}",
        "with 11881 regions",
    ]

    # When
    subprocess_result = subprocess.run(
        cmd, shell=True, capture_output=True, check=False, universal_newlines=True
    )

    # Then
    assert subprocess_result.returncode == 0
    for line in expected_stderr_substrings:
        assert line in subprocess_result.stderr
    assert output_file.exists()
    actual_first_line = output_file.read_text().strip().split("\n")[0]
    assert actual_first_line == expected_first_line
