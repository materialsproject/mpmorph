from pymatgen.core.structure import Structure

from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs
from .core import M3GNetMDMaker
from .pv_from_calc import PVFromM3GNet
from .equilibrate_volume import EquilibriumVolumeSearchMaker


def get_volume_at_temp_m3gnet_job(structure: Structure, md_inputs: M3GNetMDInputs):
    pv_md_maker = PVFromM3GNet(parameters=md_inputs)
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)
    return eq_vol_maker.make(structure)
