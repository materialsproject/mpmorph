from atomate2.vasp.jobs.core import MDMaker
from jobflow import Flow, Maker
from mpmorph.jobs.core import M3GNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromM3GNet, PVFromVasp
from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs

def get_md_flow_m3gnet(structure, temp, steps, converge_first = True, initial_vol_scale = 1):
    inputs = M3GNetMDInputs(
        temperature=temp,
        steps=steps
    )

    m3gnet_maker = M3GNetMDMaker(parameters = inputs)
    pv_md_maker = PVFromM3GNet(parameters=inputs)
    return get_md_flow(
        pv_md_maker=pv_md_maker,
        production_md_maker=m3gnet_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale
    )

def get_md_flow_vasp(structure, temperature, steps, converge_first = True, initial_vol_scale = 1):
    # TODO: Fix all of this (e.g. steps default)
    production_vasp_maker = MDMaker()
    pv_md_maker = PVFromVasp(maker=production_vasp_maker)
    production_vasp_maker.input_set_generator.user_incar_settings["NSW"] = steps

    return get_md_flow(
        pv_md_maker=pv_md_maker,
        production_md_maker=production_vasp_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale
    )


def get_md_flow(pv_md_maker, production_md_maker, structure, converge_first, initial_vol_scale):

    struct = structure.copy()
    if initial_vol_scale is not None:
        struct.scale_lattice(struct.volume * initial_vol_scale)

    if converge_first:
        return get_converge_flow(struct, pv_md_maker, production_md_maker)
    else:
        return Flow([production_md_maker.make(struct)])

def get_converge_flow(structure: Structure, pv_md_maker: Maker, production_run_maker: Maker):
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)

    equil_vol_job = eq_vol_maker.make(structure)

    final_md_job = production_run_maker.make(equil_vol_job.output)

    flow = Flow([equil_vol_job, final_md_job], output=final_md_job.output)

    return flow
