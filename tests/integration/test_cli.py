import subprocess
import shlex
import shutil

import pytest

import fur_cnvkit
from fur_cnvkit import constants

MODULE_NAME = fur_cnvkit.__name__
PROGRAM_NAME = constants.PROGRAM_NAME

# HELPERS


def get_subprocess_message(subproces_result: subprocess.CompletedProcess) -> str:
    indent = " " * 2

    msg = (
        f"Error running CLI command. "
        f"{indent}Command: {subproces_result.args}\n"
        f"{indent}Return code: {subproces_result.returncode}\n"
        f"{indent}Stdout: {subproces_result.stdout!r}\n"
        f"{indent}Stderr: {subproces_result.stderr!r}"
    )
    return msg


# TESTS


def test_python_dash_m__version():
    # Precondition
    # We assume the system has python installed but occasionally the binary may
    # be named python3 with no python binary.
    python_exec = "python"
    if shutil.which("python") is None:
        assert (
            shutil.which("python3") is not None
        ), "Python is not installed or not in PATH"
        python_exec = "python3"

    # Given
    cmd = f"{python_exec} -m {MODULE_NAME} --version"
    expected_version = fur_cnvkit.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


def test_cli_on_path():
    # When
    result = shutil.which(PROGRAM_NAME)

    assert result is not None, f"{PROGRAM_NAME} is not in PATH, has the name changed?"


def test_cli__version():
    # Given
    cmd = f"{PROGRAM_NAME} --version"
    expected_version = fur_cnvkit.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


def test_cli__help():
    # Given
    cmd = f"{PROGRAM_NAME} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout


@pytest.mark.parametrize(
    "command",
    [
        pytest.param(
            constants.COMMAND_NAME__CALCULATE_MAD,
            id="calculate_mad",
        ),
        pytest.param(
            constants.COMMAND_NAME__GENERATE_STATIC_FILES,
            id="generate_static_files",
        ),
        pytest.param(
            constants.COMMAND_NAME__GENERATE_CN_REFERENCE,
            id="generate_cn_reference",
        ),
        pytest.param(
            constants.COMMAND_NAME__GENERATE_ONCOPRINT,
            id="generate_oncoprint",
        ),
        pytest.param(
            constants.COMMAND_NAME__RUN_CNVKIT_CN_CALLING_PIPELINE,
            id="run_cnvkit_cn_calling_pipeline",
        ),
    ],
)
def test_cli__command__help(command: str):
    # Given
    cmd = f"{PROGRAM_NAME} {command} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert command in subproces_result.stdout
