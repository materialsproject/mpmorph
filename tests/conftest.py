import glob
from pathlib import Path

import pytest

TEST_FILES_PATH = Path(__file__).parent / "test_files"


@pytest.fixture(scope="session")
def md_run():
    p = glob.glob(str(TEST_FILES_PATH / "liquid_Na/Na_df*/run0/XDATCAR.gz"))
    return p
