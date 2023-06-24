from atomate2.vasp.jobs.core import MDMaker
from atomate2.vasp.sets.core import MDSetGenerator
from pymatgen.core.structure import Structure
from mpmorph.jobs.pv_from_calc import PVFromVasp

from .utils import get_md_flow
from pymatgen.io.vasp.inputs import Kpoints

VASP_MD_CONVERGE_FLOW = "VASP_MD_CONVERGE_FLOW"

def get_md_flow_vasp(
    structure: Structure,
    temperature: int,
    steps_prod: int,
    steps_pv: int,
    converge_first: bool = True,
    initial_vol_scale: int = 1
):
    my_kpoints = Kpoints()

    flow_name = f'MD_FLOW_{structure.composition.to_pretty_string()}'
    production_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=steps_prod,
            time_step=2,
            user_incar_settings={
                "ISPIN": 1, # Do not consider magnetism in AIMD simulations
                "LREAL": True, # Peform calculation in real space for AIMD due to large unit cell size
                "LAECHG": False, # Don't need AECCAR for AIMD
                "EDIFFG": None, # Does not apply to MD simulations, see: https://www.vasp.at/wiki/index.php/EDIFFG
                "GGA": "PS", # Just let VASP decide based on POTCAR - the default, PS yields the error below
                "LPLANE": False # LPLANE is recommended to be False on Cray machines (https://www.vasp.at/wiki/index.php/LPLANE)
            },
            user_kpoints_settings=my_kpoints
        ),
        # This is recommended by VASP, but defaults to false in atomate2...
        # run_vasp_kwargs={
        #     "vasp_job_kwargs": {
        #         "auto_npar": True
        #     }
        # }
    )

    pv_vasp_maker = MDMaker(
        input_set_generator=MDSetGenerator(
            ensemble="nvt",
            start_temp=temperature,
            end_temp=temperature,
            nsteps=steps_pv,
            time_step=2,
            user_incar_settings={
                "ISPIN": 1, # Do not consider magnetism in AIMD simulations
                "LREAL": True, # Peform calculation in real space for AIMD due to large unit cell size
                "LAECHG": False, # Don't need AECCAR for AIMD
                "EDIFFG": None, # Does not apply to MD simulations, see: https://www.vasp.at/wiki/index.php/EDIFFG
                "GGA": "PS", # Just let VASP decide based on POTCAR - the default, PS yields the error below
                "LPLANE": False # LPLANE is recommended to be False on Cray machines (https://www.vasp.at/wiki/index.php/LPLANE)
            },
            user_kpoints_settings=my_kpoints
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
