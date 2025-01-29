import pytest
from pathlib import Path
import json
from tempfile import TemporaryDirectory

import pandas as pd

from fur_cnvkit.generate_cnvkit_static_files import generate_parameter_file


def test_generate_parameter_file():
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
        access_bed = tempdir_path / "access.bed"
        unplaced_contig_prefixes = ["chrUn", "chrAlt"]
        parameter_file_name = "test"
        outdir = tempdir_path

        # Create mock files
        for file in [
            reference_fasta,
            baitset_bed,
            refflat_file,
            sample_metadata_xlsx,
            targets_bed,
            antitargets_bed,
            access_bed,
        ]:
            file.touch()

        # Create mock BAM files and sample metadata
        tumour_bam = tempdir_path / "tumour_sample.bam"
        normal_bam = tempdir_path / "normal_sample.bam"
        tumour_bam.touch()
        normal_bam.touch()
        bam_files = [tumour_bam, normal_bam]

        # Write mock sample metadata
        sample_metadata_content = {
            "Sanger DNA ID": ["tumour_sample", "normal_sample"],
            "T/N": ["T", "N"],
        }
        sample_metadata_df = pd.DataFrame(sample_metadata_content)
        sample_metadata_df.to_excel(sample_metadata_xlsx, index=False)

        # Call the function
        output_parameter_file = generate_parameter_file(
            bam_files=bam_files,
            reference_fasta=reference_fasta,
            unplaced_contig_prefixes=unplaced_contig_prefixes,
            baitset_bed=baitset_bed,
            refflat_file=refflat_file,
            sample_metadata_xlsx=sample_metadata_xlsx,
            access_bed=access_bed,
            targets_bed=targets_bed,
            antitargets_bed=antitargets_bed,
            parameter_file_name=parameter_file_name,
            outdir=outdir,
        )

        # Verify the output file path
        expected_output_parameter_file = (
            outdir / f"{parameter_file_name}.parameters.json"
        )
        assert output_parameter_file == expected_output_parameter_file
        assert output_parameter_file.exists()

        # Verify the contents of the generated parameter file
        with output_parameter_file.open() as json_file:
            parameter_file_data = json.load(json_file)

        expected_parameter_file_data = {
            "all_bams": [str(tumour_bam), str(normal_bam)],
            "tumour_bams": [str(tumour_bam)],
            "normal_bams": [str(normal_bam)],
            "reference_fasta": str(reference_fasta),
            "unplaced_contig_prefixes": unplaced_contig_prefixes,
            "baitset_bed": str(baitset_bed),
            "refflat_file": str(refflat_file),
            "sample_metadata_xlsx": str(sample_metadata_xlsx),
            "access_bed": str(access_bed),
            "targets_bed": str(targets_bed),
            "antitargets_bed": str(antitargets_bed),
        }

        assert parameter_file_data == expected_parameter_file_data
