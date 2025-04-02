from importlib import resources
from pathlib import Path


def get_example_sample_metadata_xlsx() -> Path:
    """
    Return the path to a xlsx file containing example sample metadata
    """
    file = resources.files("tests.mocks").joinpath("example_sample_metadata.xlsx")
    as_path = Path(str(file))
    if not as_path.exists():
        raise FileNotFoundError(f"File not found: {as_path}")
    return as_path
