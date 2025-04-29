import argparse
import sys
import fur_cnvkit
from fur_cnvkit import (
    calculate_mad,
    generate_cnvkit_static_files,
    generate_copy_number_reference,
    generate_oncoprint,
    run_cnvkit_copy_number_calling_pipeline,
)
from fur_cnvkit import constants
from fur_cnvkit.utils import logging_utils


def main():
    # Top-level parser
    parser = argparse.ArgumentParser(
        prog=constants.PROGRAM_NAME,
        description=constants.DESCRIPTION__PROGRAM,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {fur_cnvkit.__version__}",
    )

    # Command initialization
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    # Parsers for each command are set up by their respective modules
    _ = calculate_mad.get_argparser(subparser=subparsers)
    _ = generate_cnvkit_static_files.get_argparser(subparser=subparsers)
    _ = generate_copy_number_reference.get_argparser(subparser=subparsers)
    _ = generate_oncoprint.get_argparser(subparser=subparsers)
    _ = run_cnvkit_copy_number_calling_pipeline.get_argparser(subparser=subparsers)

    # Parse the arguments
    args = parser.parse_args()

    logging_utils.setup_logging(
        level=getattr(args, "verbose", "INFO"),
    )

    # Call the appropriate function based on the command from their respective
    # modules
    match args.command:
        case calculate_mad.COMMAND_NAME:
            exit_code = calculate_mad.main(args)
        case generate_cnvkit_static_files.COMMAND_NAME:
            exit_code = generate_cnvkit_static_files.main(args)
        case generate_copy_number_reference.COMMAND_NAME:
            exit_code = generate_copy_number_reference.main(args)
        case generate_oncoprint.COMMAND_NAME:
            exit_code = generate_oncoprint.main(args)
        case run_cnvkit_copy_number_calling_pipeline.COMMAND_NAME:
            exit_code = run_cnvkit_copy_number_calling_pipeline.main(args)
        case _:
            # Unlikley to see this error as the subparsers as
            # parser.parse_args() will catch unknown commands
            raise argparse.ArgumentError(
                argument=None,
                message=f"Invalid command - {args.command}",
            )

    # Exit with the exit code from the called function
    if exit_code is not None:
        sys.exit(exit_code)

    return None
