from atomate2.vasp.jobs.core import MDMaker
from atomate2.vasp.sets.core import MDSetGenerator
from pymatgen.core.structure import Structure
from mpmorph.jobs.pv_from_calc import PVFromVasp

from .utils import get_md_flow

VASP_MD_CONVERGE_FLOW = "VASP_MD_CONVERGE_FLOW"

def get_md_flow_vasp(
    structure: Structure,
    temperature: int,
    steps_prod: int,
    steps_pv: int,
    converge_first: bool = True,
    initial_vol_scale: int = 1
):
    flow_name = f'MD_FLOW_{structure.composition.to_pretty_string()}'
    production_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=steps_prod,
            time_step=2,
            user_incar_settings={
                "ISPIN": 1
            }
        )
    )

    pv_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=steps_pv,
            time_step=2,
            user_incar_settings={
                "ISPIN": 1
            }
        )
    )

    pv_extractor = PVFromVasp()

    return get_md_flow(
        pv_md_maker=pv_vasp_maker,
        pv_extractor=pv_extractor,
        production_md_maker=production_vasp_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale,
        flow_name=flow_name
    )
