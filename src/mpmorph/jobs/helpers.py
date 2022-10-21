from pymatgen.core.structure import Structure
from .core import M3GNetMDMaker
from .pv_from_calc import PVFromM3GNet
from .equilibrate_volume import EquilibriumVolumeSearchMaker

def get_volume_at_temp_m3gnet_job(structure: Structure,
                                    temp: int,
                                    steps: int = 1000):
    md_maker = M3GNetMDMaker(
        temperature = temp,
        steps = steps
    )
    pv_md_maker = PVFromM3GNet(
        md_maker = md_maker
    )                          
    eq_vol_maker = EquilibriumVolumeSearchMaker(
        pv_md_maker = pv_md_maker
    )
    
    vol_job = eq_vol_maker.make(structure)
    
    return vol_job