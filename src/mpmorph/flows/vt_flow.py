from jobflow import Flow, job
import json
import uuid

from .md_flow import get_equil_vol_flow, get_equil_vol_flow_lammps
from ..jobs.pv_from_calc import m3gnet_calc_to_vol
import pandas as pd

VOLUME_TEMPERATURE_SWEEP = "VOLUME_TEMPERATURE_SWEEP"


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
        job = get_equil_vol_flow(structure=structure, temp=temp, steps=steps)
        volume_jobs.append(job)
        vs.append(job.output.volume)

    collect_job = _collect_vt_results(vs, temps, structure, output_name)

    new_flow = Flow(
        [*volume_jobs, collect_job],
        output=collect_job.output,
        name=VOLUME_TEMPERATURE_SWEEP,
    )
    return new_flow


def get_vt_sweep_flow_lammps(
    structure,
    lower_bound=100,
    upper_bound=1100,
    temp_step=100,
    output_name="vt.out",
    steps=2000,
    mp_id=None,
):
    v_outputs = []
    volume_jobs = []
    temps = list(range(lower_bound, upper_bound, temp_step))

    for temp in temps:
        job = get_equil_vol_flow_lammps(
            structure=structure,
            temp=temp,
            steps=steps,
        )
        volume_jobs.append(job)
        v_outputs.append(job.output.output)

    collect_job = _collect_vt_results(v_outputs, temps, structure, output_name, mp_id)

    flow_name = f"{structure.composition.reduced_formula}-Melting Point"
    new_flow = Flow(
        [*volume_jobs, collect_job], output=collect_job.output, name=flow_name
    )
    return new_flow


@job
def _collect_vt_results(v_outputs, ts, structure, output_fn, mp_id):
    result = {
        "structure": structure.as_dict(),
        "volumes": [get_converged_vol(v) for v in v_outputs],
        "temps": ts,
        "mp_id": mp_id,
        "reduced_formula": structure.composition.reduced_formula,
        "formula": structure.composition.formula,
        "uuid": str(uuid.uuid4()),
    }

    with open(output_fn, "+w") as f:
        f.write(json.dumps(result))
    return result


def get_converged_vol(v_output):
    df = pd.DataFrame.from_dict(v_output)
    total_steps = (len(df) - 1) * 10
    avging_window = int(total_steps / 30)
    vols = df.iloc[-avging_window::]["vol"]
    eq_vol = vols.values.mean()
    return float(eq_vol)
