import glob
from pathlib import Path

import pytest
from jobflow.core.store import JobStore
from maggma.stores import MemoryStore
from monty.serialization import loadfn

TEST_FILES_PATH = Path(__file__).parent / "test_files"
AL_STRUCTURE = loadfn(TEST_FILES_PATH / "al_structure.json.gz")


@pytest.fixture(scope="session")
def al_structure():
    return AL_STRUCTURE


@pytest.fixture(scope="session")
def md_run():
    p = glob.glob(str(TEST_FILES_PATH / "liquid_Na/Na_df*/run0/XDATCAR.gz"))
    return p


@pytest.fixture(scope="session")
def job_store():
    additional_stores = {"trajectory": MemoryStore()}
    return JobStore(MemoryStore(), additional_stores=additional_stores)
