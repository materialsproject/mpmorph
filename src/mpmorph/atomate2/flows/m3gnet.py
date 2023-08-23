from jobflow import Flow
from ..jobs.core import M3GNetMDMaker

from ..jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from ..jobs.pv_from_calc import PVFromM3GNet
from ..jobs.tasks.m3gnet_input import M3GNetMDInputs
from .utils import get_md_flow, collect_vt_results, VOLUME_TEMPERATURE_SWEEP

EQUILIBRATE_VOLUME_FLOW = "EQUILIBRATE_VOLUME_FLOW"
M3GNET_MD_FLOW = "M3GNET_MD_FLOW"
M3GNET_MD_CONVERGED_VOL_FLOW = "M3GNET_MD_CONVERGED_VOL_FLOW"
LAMMPS_VOL_FLOW = "LAMMPS_VOL_FLOW"

def get_md_flow_m3gnet(
        structure: Structure,
        temp: int,
        steps_prod: int,
        steps_pv: int,
        converge_first: bool = True,
        initial_vol_scale: float = 1,
        **input_kwargs
    ):
    inputs_prod = M3GNetMDInputs(temperature=temp, steps=steps_prod, **input_kwargs)
    inputs_pv = M3GNetMDInputs(temperature=temp, steps=steps_pv, **input_kwargs)

    pv_extractor = PVFromM3GNet()
    
    pv_m3gnet_maker = M3GNetMDMaker(parameters=inputs_pv)
    prod_m3gnet_maker = M3GNetMDMaker(parameters=inputs_prod)

    return get_md_flow(
        pv_md_maker=pv_m3gnet_maker,
        pv_extractor=pv_extractor,
        production_md_maker=prod_m3gnet_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale,
    )


def get_equil_vol_flow(
        structure: Structure,
        temp: int,
        steps: int
    ):
    inputs = M3GNetMDInputs(temperature=temp, steps=steps)

    pv_extractor = PVFromM3GNet()
    pv_md_maker = M3GNetMDMaker(parameters=inputs)

    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker, pv_extractor=pv_extractor)
    equil_vol_job = eq_vol_maker.make(structure)
    flow = Flow(
        [equil_vol_job], output=equil_vol_job.output, name=EQUILIBRATE_VOLUME_FLOW
    )
    return flow

def get_vt_sweep_flow(
    structure: Structure,
    lower_bound: int = 100,
    upper_bound: int = 1100,
    temp_step: int = 100,
    output_name: str = "vt.out",
    steps: int = 2000,
):
    vs = []
    volume_jobs = []
    temps = list(range(lower_bound, upper_bound, temp_step))

    for temp in temps:
        job = get_equil_vol_flow(structure=structure, temp=temp, steps=steps)
        volume_jobs.append(job)
        vs.append(job.output.volume)

    collect_job = collect_vt_results(vs, temps, structure, output_name)

    new_flow = Flow(
        [*volume_jobs, collect_job],
        output=collect_job.output,
        name=VOLUME_TEMPERATURE_SWEEP,
    )
    return new_flow