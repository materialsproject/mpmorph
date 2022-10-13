import pytest
from jobflow.managers.local import run_locally

from mpmorph.jobs.core import M3GNetMDMaker


@pytest.fixture
def m3gnet_job(al_structure):
    return M3GNetMDMaker().make(al_structure, steps=10)


def test_m3gnet_job(m3gnet_job, job_store):
    output = run_locally(m3gnet_job, store=job_store, ensure_success=True)

    d = output[m3gnet_job.uuid][1].output
    assert d.__class__.__name__ == "M3GNetMDCalculation"
