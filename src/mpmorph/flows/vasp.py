from atomate2.vasp.jobs.core import MDMaker
from atomate2.vasp.sets.core import MDSetGenerator
from pymatgen.core.structure import Structure
from mpmorph.jobs.pv_from_calc import PVFromVasp

from .utils import get_md_flow

VASP_MD_CONVERGE_FLOW = "VASP_MD_CONVERGE_FLOW"

def get_md_flow_vasp(
    structure: Structure,
    temperature: int,
    production_steps: int,
    equilibration_steps: int,
    converge_first=True,
    initial_vol_scale=1,
):
    production_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=production_steps,
            time_step=2,
        )
    )

    pv_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=equilibration_steps,
            time_step=2,
        )
    )

    pv_md_maker = PVFromVasp(md_maker=pv_vasp_maker)

    return get_md_flow(
        pv_md_maker=pv_md_maker,
        production_md_maker=production_vasp_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale,
        flow_name=VASP_MD_CONVERGE_FLOW
    )


