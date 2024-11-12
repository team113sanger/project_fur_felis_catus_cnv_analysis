import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from utils.cnvkit_utils import run_command, run_cnvkit_access, run_cnvkit_autobin


# Tests for run_cnvkit_access
@patch("utils.cnvkit_utils.run_command")
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
@patch("utils.cnvkit_utils.run_command")
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
