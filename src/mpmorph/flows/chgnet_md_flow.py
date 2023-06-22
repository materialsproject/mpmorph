from jobflow import Flow
from mpmorph.jobs.core import CHGNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromCHGNet
from mpmorph.jobs.tasks.chgnet_input import CHGNetMDInputs
from .utils import get_md_flow

EQUILIBRATE_VOLUME_FLOW = "EQUILIBRATE_VOLUME_FLOW"
CHGNET_MD = "CHGNET_MD"

def get_md_flow_chgnet(
        structure: Structure,
        temp: int,
        steps_prod: int,
        steps_pv: int,
        converge_first: bool = True,
        initial_vol_scale: float = 1,
        **input_kwargs
    ):
    inputs_prod = CHGNetMDInputs(temperature=temp, steps=steps_prod, **input_kwargs)
    inputs_pv = CHGNetMDInputs(temperature=temp, steps=steps_pv, **input_kwargs)

    pv_extractor = PVFromCHGNet()
    
    pv_m3gnet_maker = CHGNetMDMaker(parameters=inputs_pv)
    prod_m3gnet_maker = CHGNetMDMaker(parameters=inputs_prod)

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
    inputs = CHGNetMDInputs(temperature=temp, steps=steps)

    pv_extractor = PVFromCHGNet()
    pv_md_maker = CHGNetMDMaker(parameters=inputs)

    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker, pv_extractor=pv_extractor)
    equil_vol_job = eq_vol_maker.make(structure)
    flow = Flow(
        [equil_vol_job], output=equil_vol_job.output, name=EQUILIBRATE_VOLUME_FLOW
    )
    return flow
