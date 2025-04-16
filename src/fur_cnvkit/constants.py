"""
This module contains string constants used in the CLI, centralizing
them in one place for easy maintenance.
"""

# Program name
#
# The name of the program that appears in the unified CLI entrypoint help message.
PROGRAM_NAME: str = "fur_cnvkit"

# Pacakged R script
PACKAGED_R_SCRIPT: str = "generate_copy_number_heatmap_plot.R"

# Command names
#
# Command names appear in the help message and are used to determine which
# subcommand to run.
COMMAND_NAME__CALCULATE_MAD: str = "mad"
COMMAND_NAME__GENERATE_STATIC_FILES: str = "cnvkit_static_files"
COMMAND_NAME__GENERATE_CN_REFERENCE: str = "copy_number_reference"
COMMAND_NAME__GENERATE_ONCOPRINT: str = "oncoprint"
COMMAND_NAME__RUN_CNVKIT_CN_CALLING_PIPELINE: str = "cnvkit_cn_calling_pipeline"

# Program & Command descriptions
#
# Text that appears in the help message that describes a command or the program
# itself.
DESCRIPTION__PROGRAM: str = "FUR CNVkit - a collection commands for CN analysis as part of FUR project that rely on CNVkit."
DESCRIPTION__CALCULATE_MAD: str = "Run the weighted Median Absolute Deviation (MAD) vs. segment size calculation pipeline on CNVkit files."
DESCRIPTION__GENERATE_STATIC_FILES: str = (
    "Generate various CNVKit files needed to run downstream analyses."
)
DESCRIPTION__GENERATE_CN_REFERENCE: str = (
    "Generate a copy number reference file for a given cohort of samples."
)
DESCRIPTION__GENERATE_ONCOPRINT: str = "Generate an oncoprint figure showing recurrent somatic mutations and CNVs across samples."
DESCRIPTION__RUN_CNVKIT_CN_CALLING_PIPELINE: str = (
    "Run the end-to-end copy number calling workflow for a cohort of samples."
)

# Command short-help
#
# Text that appears in the program help that concisely describes the command.
SHORT_HELP__CALCULATE_MAD: str = (
    "Calculate the weighted Median Absolute Deviation (MAD) of segments in samples."
)
SHORT_HELP__GENERATE_STATIC_FILES: str = "Generate CNVKit static files."
SHORT_HELP__GENERATE_CN_REFERENCE: str = "Generate a copy number reference file."
SHORT_HELP__GENERATE_ONCOPRINT: str = "Generate an oncoprint figure."
SHORT_HELP__RUN_CNVKIT_CN_CALLING_PIPELINE: str = (
    "Run the end-to-end copy number calling workflow for a cohort of samples."
)
