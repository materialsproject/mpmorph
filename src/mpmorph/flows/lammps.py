from atomate2.vasp.jobs.core import MDMaker
from atomate2.vasp.sets.core import MDSetGenerator
from jobflow import Flow, Maker
from mpmorph.jobs.core import M3GNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from mpmorph.jobs.lammps.lammps_basic_const_temp import BasicLammpsConstantTempMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromCalc, PVFromM3GNet, PVFromVasp
from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs

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



