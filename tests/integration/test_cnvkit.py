import os
import pathlib
import subprocess
import re
import shutil
import pytest
import importlib.resources

from fur_cnvkit import constants
from tests import conftest


# -------------------------------
# Fixtures are defined in tests/conftest.py
# -------------------------------


# -------------------------------
# Additional helpers to determine if test data is available
# -------------------------------


def should_skip_tests() -> bool:
    """
    Returns True if required BAM and baitset files are not found.
    This can be used with @pytest.mark.skipif.
    """

    def _can_find_bams() -> bool:
        raw_bam_dir = os.environ.get(conftest.ENV_VAR_BAM_DIR, "")
        bam_dir = pathlib.Path(raw_bam_dir)
        return bool(raw_bam_dir) and bam_dir.exists() and any(bam_dir.glob("**/*.bam"))

    def _can_find_baitset() -> bool:
        raw_baitset = os.environ.get(conftest.ENV_VAR_BAITSET_DIR, "")
        baitset = pathlib.Path(raw_baitset)
        return bool(raw_baitset) and baitset.exists() and any(baitset.glob("**/*.bed"))

    return not all([_can_find_bams(), _can_find_baitset()])


# -------------------------------
# Test functions
# -------------------------------


def test_cnvkit_is_installed():
    program = "cnvkit.py"
    assert bool(shutil.which(program)), "cnvkit.py is not installed or not in PATH"


def test_cnvkit_version_is_expected():
    # Given
    expected_version = "0.9.10"
    cmd = "cnvkit.py version"

    # When
    version = subprocess.check_output(cmd, shell=True).decode().strip()

    # Then
    assert (
        version == expected_version
    ), f"Expected version {expected_version}, got {version}"


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

    cmd = (
        f"cnvkit.py autobin -f {feline_reference_fasta} -m {method} -t {feline_baitset} "
        f"--annotate {feline_refflat_file} --target-output-bed {expected_output_target_bed} "
        f"--antitarget-output-bed {expected_output_antitarget_bed} "
        f"{' '.join(str(bam) for bam in bams)}"
    )

    # When
    subprocess_result = subprocess.run(
        cmd, shell=True, capture_output=True, check=False, universal_newlines=True
    )

    # Then
    assert (
        subprocess_result.returncode == 0
    ), f"Command failed: {subprocess_result.stderr}"
    assert (
        expected_output_target_bed.exists()
    ), f"Target bed file not created: {expected_output_target_bed}"
    assert (
        expected_output_antitarget_bed.exists()
    ), f"Antitarget bed file not created: {expected_output_antitarget_bed}"
    for substring in expected_stderr_substrings:
        assert (
            substring in subprocess_result.stderr
        ), f"Expected '{substring}' in stderr"
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
    assert (
        subprocess_result.returncode == 0
    ), f"Command failed: {subprocess_result.stderr}"
    for substring in expected_stderr_substrings:
        assert (
            substring in subprocess_result.stderr
        ), f"Expected '{substring}' in stderr"
    assert output_file.exists(), f"Output file not created: {output_file}"
    actual_first_line = output_file.read_text().strip().split("\n")[0]
    assert (
        actual_first_line == expected_first_line
    ), f"Expected first line '{expected_first_line}', got '{actual_first_line}'"


def test_is_R_installed():
    # Given
    programs = ["R", "Rscript"]

    # When
    is_installed = all(bool(shutil.which(program)) for program in programs)

    # Then
    assert is_installed, "R or Rscript is not installed or not in PATH"


def test_is_our_R_script_installed():
    # Given
    import fur_cnvkit

    expected_path = (
        pathlib.Path(str(importlib.resources.files(fur_cnvkit)))
        / constants.PACKAGED_R_SCRIPT
    )
    expected_path = expected_path.resolve()

    # When
    is_extant = importlib.resources.is_resource(fur_cnvkit, constants.PACKAGED_R_SCRIPT)

    # Then
    err_msg = (
        f"The R script {constants.PACKAGED_R_SCRIPT!r} cannot be found. "
        "Has its name been changed? Or is it in a different location - expected to be found"
        f"{str(expected_path)!r}?"
    )
    assert is_extant, err_msg


@pytest.mark.parametrize(
    "package_name, used_by",
    [
        pytest.param(
            "DNAcopy",
            "cnvkit",
            id="DNAcopy",
        ),
        pytest.param(
            "ComplexHeatmap",
            "cnvkit",
            id="ComplexHeatmap",
        ),
        pytest.param(
            "ComplexHeatmap",
            f"our script ({constants.PACKAGED_R_SCRIPT})",
            id="ComplexHeatmap#2",
        ),
        pytest.param(
            "optparse",
            f"our script ({constants.PACKAGED_R_SCRIPT})",
            id="optparse",
        ),
        pytest.param(
            "dplyr",
            f"our script ({constants.PACKAGED_R_SCRIPT})",
            id="dplyr",
        ),
        pytest.param(
            "circlize",
            f"our script ({constants.PACKAGED_R_SCRIPT})",
            id="circlize",
        ),
        pytest.param(
            "randomcoloR",
            f"our script ({constants.PACKAGED_R_SCRIPT})",
            id="randomcoloR",
        ),
    ],
)
def test_R_package_is_installed(package_name: str, used_by: str):
    # Given
    from fur_cnvkit.utils import cnvkit_utils

    # When
    actually_exists = cnvkit_utils.is_R_package_installed(package_name)

    # Then
    err_msg = (
        f"The R package {package_name!r} cannot be found. "
        f"It is needed by {used_by}. "
        "Is it installed? Are the paths correct?"
    )
    assert actually_exists, err_msg
