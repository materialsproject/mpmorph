import pytest
from jobflow.managers.local import run_locally

from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs
from mpmorph.jobs.core import M3GNetMDMaker


@pytest.fixture
def m3gnet_job(al_structure):
    params = M3GNetMDInputs(
        steps = 10
    )
    print(params)
    return M3GNetMDMaker(parameters=params).make(al_structure)


def test_m3gnet_job(m3gnet_job, job_store):
    output = run_locally(m3gnet_job, store=job_store, ensure_success=True)

    d = output[m3gnet_job.uuid][1].output
    assert d.__class__.__name__ == "M3GNetMDCalculation"
