import shutil
import pathlib
import os
import subprocess
import re

import pytest

from fur_cnvkit.utils import cnvkit_utils

# CONSTANTS


ENV_VAR_BAM_DIR = "TEST_BAM_DIR"
ENV_VAR_BAITSET_DIR = "TEST_BAITSET_DIR"
ENV_VAR_GENOME_DIR = "TEST_GENOME_DIR"
ENV_VAR_ANNOTATION_DIR = "TEST_ANNOTATION_DIR"


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
def test_cnvkit_runs__numpy_regression(
    feline_reference_fasta: pathlib.Path,
    feline_baitset: pathlib.Path,
    feline_refflat_file: pathlib.Path,
    bams: list[pathlib.Path],
    tmp_path: pathlib.Path,
):
    # Given
    method = "hybrid"
    expected_output_target_bed = tmp_path / feline_baitset.name.replace(
        ".bed", ".target.bed"
    )
    expected_output_antitarget_bed = tmp_path / feline_baitset.name.replace(
        ".bed", ".antitarget.bed"
    )
    expected_stderr_substrings = [
        "Detected file format: bed",
        "Detected file format: refflat",
        f"Wrote {expected_output_target_bed}",
        f"Wrote {expected_output_antitarget_bed}",
    ]
    numpy_error_pattern = r"np\.asfarray.*was removed in the NumPy 2\.0 release"

    cmd = f"cnvkit.py autobin -f {feline_reference_fasta} -m {method} -t {feline_baitset} --annotate {feline_refflat_file} --target-output-bed {expected_output_target_bed} --antitarget-output-bed {expected_output_antitarget_bed} {' '.join([str(bam) for bam in bams])}"

    # When
    subprocess_result = subprocess.run(
        cmd, shell=True, capture_output=True, check=False, universal_newlines=True
    )

    # Then
    assert subprocess_result.returncode == 0
    assert pathlib.Path(expected_output_target_bed).exists()
    assert pathlib.Path(expected_output_antitarget_bed).exists()
    for line in expected_stderr_substrings:
        assert line in subprocess_result.stderr
    assert not re.search(
        numpy_error_pattern, subprocess_result.stderr
    ), f"Detected numpy asfarray error: {subprocess_result.stderr}"


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


def test_is_R_installed():
    # Given
    should_exist = True
    programs = ["R", "Rscript"]

    # When
    is_installed = all([bool(shutil.which(program)) for program in programs])

    # Then
    assert is_installed == should_exist


def test_is_DNACopy_a_findable_R_package():
    # Given
    package = "DNAcopy"
    should_exist = True

    # When
    actually_exists = cnvkit_utils.is_R_package_installed(package)

    # Then
    err_msg = f"The R pacakge {package!r} can not be found. Is it installed? Are the paths correct?"
    assert actually_exists == should_exist, err_msg


def test_is_ComplexHeatmap_a_findable_R_package():
    # Given
    package = "ComplexHeatmap"
    should_exist = True

    # When
    actually_exists = cnvkit_utils.is_R_package_installed(package)

    # Then
    err_msg = f"The R pacakge {package!r} can not be found. Is it installed? Are the paths correct?"
    assert actually_exists == should_exist, err_msg
