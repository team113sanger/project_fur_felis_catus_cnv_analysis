from pathlib import Path
import shutil

import pytest
from unittest.mock import patch, MagicMock

from fur_cnvkit.utils.cnvkit_utils import (
    run_command,
    run_cnvkit_access,
    run_cnvkit_autobin,
)
from fur_cnvkit.utils import cnvkit_utils


# Tests for run_cnvkit_access
@patch("fur_cnvkit.utils.cnvkit_utils.run_command")
def test_run_cnvkit_access(mock_run_command, feline_reference_fasta: Path):
    # Given
    outdir = Path("/test/outdir")
    expected_output_bed_name = Path(
        "access-Felis_catus.Felis_catus_9.0.dna.toplevel.bed"
    )
    expected_output_bed_path = outdir / expected_output_bed_name
    expected_cmd = f"cnvkit.py access {str(feline_reference_fasta)} -o {str(expected_output_bed_path)}"

    mock_run_command.return_value = None

    # When
    output_bed_path = run_cnvkit_access(feline_reference_fasta, outdir)

    # Then
    assert output_bed_path == expected_output_bed_path
    mock_run_command.assert_called_once_with(expected_cmd)


# Tests for run_cnvkit_autobin
@patch("fur_cnvkit.utils.cnvkit_utils.run_command")
def test_run_cnvkit_autobin(
    mock_run_command,
    bams: list[Path],
    feline_baitset: Path,
    feline_cnvkit_access_bed: Path,
    feline_refflat_file: Path,
):
    # Given
    outdir = Path("/test/outdir")
    expected_target_bed_name = (
        "S3250994_Feline_HSA_Jan2020_146_canonical_pad100.merged.target.bed"
    )
    expected_antitarget_bed_name = (
        "S3250994_Feline_HSA_Jan2020_146_canonical_pad100.merged.antitarget.bed"
    )
    expected_target_bed_path = outdir / expected_target_bed_name
    expected_antitarget_bed_path = outdir / expected_antitarget_bed_name
    bam_files_str = " ".join(map(str, bams))
    expected_cmd = (
        f"cnvkit.py autobin {bam_files_str} -t {str(feline_baitset)} "
        f"-g {str(feline_cnvkit_access_bed)} --annotate {feline_refflat_file} "
        f"--target-output-bed {str(expected_target_bed_path)} "
        f"--antitarget-output-bed {str(expected_antitarget_bed_path)}"
    )

    mock_run_command.return_value = None

    # When
    target_bed_dict = run_cnvkit_autobin(
        bams, feline_baitset, feline_cnvkit_access_bed, feline_refflat_file, outdir
    )

    # Then
    assert target_bed_dict["target"] == expected_target_bed_path
    assert target_bed_dict["antitarget"] == expected_antitarget_bed_path
    mock_run_command.assert_called_once_with(expected_cmd)


@pytest.mark.parametrize(
    "package, should_exist",
    [
        # Valid R stdlib-packages
        pytest.param("stats", True, id="R-stdlib-stats"),
        pytest.param("utils", True, id="R-stdlib-utils"),
        pytest.param("grDevices", True, id="R-stdlib-grDevices"),
        pytest.param("graphics", True, id="R-stdlib-graphics"),
        pytest.param("datasets", True, id="R-stdlib-datasets"),
        pytest.param("methods", True, id="R-stdlib-methods"),
        pytest.param("base", True, id="R-stdlib-base"),
        # Invalid R packages
        pytest.param("non_existent_package", False, id="invalid#non_existent_package"),
    ],
)
def test_is_R_package_installed(package: str, should_exist: bool):
    # When
    actually_exists = cnvkit_utils.is_R_package_installed(package)

    # Then
    assert actually_exists == should_exist


@pytest.mark.xfailif(
    not shutil.which("Rscript"),
    reason="Rscript not found in PATH, failing test",
)
def test_get_libPaths_used_by_R():
    # When
    libPaths = cnvkit_utils.get_libPaths_used_by_R()

    # Then
    assert isinstance(libPaths, list)
    assert all([Path(p).exists() for p in libPaths])
    assert len(libPaths) > 0, "No R library paths found, 1 expected at least"
