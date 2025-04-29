from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from fur_cnvkit.utils.cnvkit_utils import (
    run_command,
    run_cnvkit_access,
    run_cnvkit_autobin,
    perform_centring,
    filter_genemetrics_file,
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


@pytest.mark.xfail(
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


@pytest.fixture
def mock_paths():
    """Fixture to provide mock input and output paths."""
    return Path("mock_input.cns"), Path("mock_output/")


@patch("fur_cnvkit.utils.cnvkit_utils.execute_command")
def test_perform_centring_success(mock_execute_command, mock_paths):
    """Test perform_centring with a valid log2 shift output."""

    # Mock successful command output
    mock_stderr = "Shifting log2 values by 0.039127\nWrote mock_output/mock_input.median_centred.cns with 20 regions"
    mock_result = MagicMock(stderr=mock_stderr)
    mock_execute_command.return_value = mock_result

    input_file, output_dir = mock_paths
    output_file, shift_value = perform_centring(input_file, output_dir)

    assert output_file == output_dir / "mock_input.median_centred.cns"
    assert shift_value == 0.039127  # Extracted from the mocked output


@patch("fur_cnvkit.utils.cnvkit_utils.execute_command")
def test_perform_centring_no_shift_value(mock_execute_command, mock_paths):
    """Test perform_centring when no shift value is found."""

    # Mock output without log2 shift value
    mock_stderr = "Wrote mock_output/mock_input.median_centred.cns with 20 regions"
    mock_result = MagicMock(stderr=mock_stderr)
    mock_execute_command.return_value = mock_result

    input_file, output_dir = mock_paths
    output_file, shift_value = perform_centring(input_file, output_dir)

    assert output_file == output_dir / "mock_input.median_centred.cns"
    assert shift_value is None  # No shift value in output


@patch("fur_cnvkit.utils.cnvkit_utils.execute_command")
def test_perform_centring_command_failure(mock_execute_command, mock_paths):
    """Test perform_centring when the command fails."""

    # Mock subprocess failure
    mock_execute_command.side_effect = subprocess.CalledProcessError(
        1, "cnvkit.py call"
    )

    input_file, output_dir = mock_paths
    output_file, shift_value = perform_centring(input_file, output_dir)

    assert output_file == output_dir / "mock_input.median_centred.cns"
    assert shift_value is None  # Expect None due to command failure


@pytest.fixture
def genemetrics_test_file():
    """Fixture to create a temporary genemetrics test file."""
    temp_dir = Path(tempfile.mkdtemp())
    test_file = temp_dir / "test.genemetrics.out"

    # Sample genemetrics data
    data = """gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tprobes
    GeneA\tchr1\t100\t200\t-1.5\t30\t0.8\t50
    GeneB\tchr2\t300\t400\t0.5\t20\t0.7\t40
    GeneC\tchr3\t500\t600\t2.0\t25\t0.6\t30
    GeneD\tchr4\t700\t800\t-3.0\t35\t0.9\t60
    """

    # Write test data to file
    with open(test_file, "w") as f:
        f.write(data)

    return test_file, temp_dir


def test_filter_genemetrics_file(genemetrics_test_file):
    """Test filtering function with sample data."""
    genemetrics_file, output_directory = genemetrics_test_file
    lower_threshold = -2.0
    upper_threshold = 1.5

    # Call the function
    output_file = filter_genemetrics_file(
        genemetrics_file, lower_threshold, upper_threshold, output_directory
    )

    # Verify output file exists
    assert output_file.exists(), "Filtered output file was not created"

    # Read the output file
    df_filtered = pd.read_csv(output_file, sep="\t")

    # Strip any leading/trailing spaces from 'gene' column
    df_filtered["gene"] = df_filtered["gene"].str.strip()

    # Expected filtered results (only GeneA, GeneC, and GeneD should be retained)
    expected_genes = {"GeneC", "GeneD"}

    print(df_filtered)

    # Validate expected genes are in the output
    assert (
        set(df_filtered["gene"]) == expected_genes
    ), "Filtered results do not match expected genes"
