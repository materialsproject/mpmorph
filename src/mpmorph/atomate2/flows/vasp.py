from atomate2.vasp.jobs.md import MDMaker
from atomate2.vasp.sets.core import MDSetGenerator
from pymatgen.core.structure import Structure
from ..jobs.pv_from_calc import PVFromVasp


from .utils import get_md_flow
from pymatgen.io.vasp.inputs import Kpoints

from typing import Dict

VASP_MD_CONVERGE_FLOW = "VASP_MD_CONVERGE_FLOW"

def get_md_flow_vasp(
    structure: Structure,
    converge_first: bool = True,
    temperature: int = None,
    steps_prod: int = None,
    steps_pv: int = None,
    initial_vol_scale: int = 1,
    scale_factor_increment: float = 0.2,
    production_md_user_incar_settings: Dict = None,
    pv_user_incar_settings: Dict = None
):

    incar_settings = {
        "ISPIN": 1, # Do not consider magnetism in AIMD simulations
        "LREAL": "Auto", # Peform calculation in real space for AIMD due to large unit cell size
        "LAECHG": False, # Don't need AECCAR for AIMD
        "EDIFFG": None, # Does not apply to MD simulations, see: https://www.vasp.at/wiki/index.php/EDIFFG
        "GGA": "PS", # Just let VASP decide based on POTCAR - the default, PS yields the error below
        "LPLANE": False, # LPLANE is recommended to be False on Cray machines (https://www.vasp.at/wiki/index.php/LPLANE)
        "LDAUPRINT": 0,
    }

    gamma_point = Kpoints()

    if production_md_user_incar_settings is None:
        production_md_user_incar_settings = {}

    production_md_user_incar_settings = { **incar_settings, **production_md_user_incar_settings }

    production_md_set_generator = MDSetGenerator(
        ensemble="nvt",
        start_temp=temperature,
        end_temp=temperature,
        nsteps=steps_prod,
        time_step=2,
        user_incar_settings=production_md_user_incar_settings,
        user_kpoints_settings=gamma_point
    )

    if pv_user_incar_settings is None:
        pv_user_incar_settings = {}

    pv_user_incar_settings = { **incar_settings, **pv_user_incar_settings }

    pv_md_set_generator = MDSetGenerator(
        ensemble="nvt",
        start_temp=temperature,
        end_temp=temperature,
        nsteps=steps_pv,
        time_step=2,
        user_incar_settings=pv_user_incar_settings,
        user_kpoints_settings=gamma_point
    )


    flow_name = f'MD_FLOW_{structure.composition.to_pretty_string()}'
    production_vasp_maker = MDMaker(
        input_set_generator=production_md_set_generator
    )

    pv_vasp_maker = MDMaker(
        input_set_generator=pv_md_set_generator
    )

    pv_extractor = PVFromVasp()

    return get_md_flow(
        pv_md_maker=pv_vasp_maker,
        pv_extractor=pv_extractor,
        production_md_maker=production_vasp_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale,
        scale_factor_increment=scale_factor_increment,
        flow_name=flow_name
    )
