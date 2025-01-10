import pytest
from pathlib import Path
import json
from tempfile import TemporaryDirectory
from fur_cnvkit.generate_cnvkit_static_files import generate_config_file


def test_generate_config_file():
    # Create temporary files and directories for the test
    with TemporaryDirectory() as tempdir:
        tempdir_path = Path(tempdir)

        # Set up test inputs
        reference_fasta = tempdir_path / "reference.fasta"
        baitset_bed = tempdir_path / "baitset.bed"
        refflat_file = tempdir_path / "refflat.txt"
        sample_metadata_xlsx = tempdir_path / "sample_metadata.xlsx"
        targets_bed = tempdir_path / "targets.bed"
        antitargets_bed = tempdir_path / "antitargets.bed"
        config_file_name = "test"
        outdir = tempdir_path

        # Create mock files
        for file in [
            reference_fasta,
            baitset_bed,
            refflat_file,
            targets_bed,
            antitargets_bed,
        ]:
            file.touch()

        # Call the function
        output_config = generate_config_file(
            reference_fasta=reference_fasta,
            baitset_bed=baitset_bed,
            refflat_file=refflat_file,
            sample_metadata_xlsx=sample_metadata_xlsx,
            targets_bed=targets_bed,
            antitargets_bed=antitargets_bed,
            config_file_name=config_file_name,
            outdir=outdir,
        )

        # Verify the output file path
        expected_output_config = outdir / f"{config_file_name}.config.json"
        assert output_config == expected_output_config
        assert output_config.exists()

        # Verify the contents of the generated config file
        with output_config.open() as json_file:
            config_data = json.load(json_file)

        expected_config_data = {
            "reference_fasta": str(reference_fasta),
            "baitset_bed": str(baitset_bed),
            "refflat_file": str(refflat_file),
            "sample_metadata_xlsx": str(sample_metadata_xlsx),
            "targets_bed": str(targets_bed),
            "antitargets_bed": str(antitargets_bed),
        }

        assert config_data == expected_config_data
