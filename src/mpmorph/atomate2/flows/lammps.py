from jobflow import Flow

from mpmorph.jobs.lammps.lammps_basic_const_temp import BasicLammpsConstantTempMaker
from pymatgen.core.structure import Structure

from .utils import collect_vt_results

LAMMPS_VOL_FLOW = "LAMMPS_VOL_FLOW"


def get_equil_vol_flow_lammps(structure: Structure,
                              temp: int,
                              steps: int):
    vol_maker = BasicLammpsConstantTempMaker()
    vol_job = vol_maker.make(
        temp,
        steps,
        structure
    )
    flow = Flow([vol_job], output=vol_job, name=LAMMPS_VOL_FLOW)
    return flow


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

    collect_job = collect_vt_results(v_outputs, temps, structure, output_name, mp_id)

    flow_name = f"{structure.composition.reduced_formula}-Melting Point"
    new_flow = Flow(
        [*volume_jobs, collect_job], output=collect_job.output, name=flow_name
    )
    return new_flow