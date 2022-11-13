from jobflow import Flow, job
import json

from .md_flow import get_md_flow_m3gnet
from ..jobs.pv_from_calc import m3gnet_calc_to_vol

def get_vt_sweep_flow(
    structure,
    lower_bound=100,
    upper_bound=1100,
    temp_step=100,
    output_name="vt.out",
    steps=2000,
):

    vs = []
    volume_jobs = []
    temps = list(range(lower_bound, upper_bound, temp_step))

    for temp in temps:
        job = get_md_flow_m3gnet(
            structure=structure,
            temp=temp,
            steps=steps,
            converge_first=True
        )
        volume_jobs.append(job)
        vs.append(m3gnet_calc_to_vol(job.output))

    collect_job = _collect_vt_results(vs, temps, structure, output_name)

    new_flow = Flow([*volume_jobs, collect_job], output=collect_job.output)
    return new_flow


@job
def _collect_vt_results(vs, ts, structure, output_fn):
    result = {"structure": structure.as_dict(), "volumes": vs, "temps": ts}

    with open(output_fn, "+w") as f:
        f.write(json.dumps(result))
    return result
