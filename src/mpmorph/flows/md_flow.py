from atomate2.vasp.jobs.core import MDMaker
from jobflow import Flow, Maker
from mpmorph.jobs.core import M3GNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromM3GNet, PVFromVasp

# TODO: Fix all of this (e.g. steps default)
def get_converge_flow_vasp(structure, temperature, steps=10):
    md_maker = MDMaker()
    pv_md_maker = PVFromVasp(maker=md_maker)
    md_maker.input_set_generator.user_incar_settings["NSW"] = steps
    return get_converge_flow(structure, pv_md_maker)


def get_converge_flow_m3gnet(structure, temp, steps: int = 1000):
    md_maker = M3GNetMDMaker(temperature=temp, steps=steps)
    pv_md_maker = PVFromM3GNet(md_maker=md_maker)
    return get_converge_flow(structure, pv_md_maker)


def get_converge_flow(structure: Structure, pv_md_maker: Maker):
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)

    equil_vol_job = eq_vol_maker.make(structure)

    final_md_job = pv_md_maker.md_maker.make(equil_vol_job.output)

    flow = Flow([equil_vol_job, final_md_job], output=final_md_job.output)

    return flow
