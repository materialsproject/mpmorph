
import pytest
from pymatgen.io.vasp.outputs import Xdatcar

from mpmorph.analysis.diffusion import Diffusion


@pytest.fixture
def diffusion(md_run):
    structures = Xdatcar(md_run[2]).structures
    d = Diffusion(structures, 300, 2, skip_first=250)
    return d


def test_diffusion(diffusion):
    assert diffusion.getD("Na") is not None
