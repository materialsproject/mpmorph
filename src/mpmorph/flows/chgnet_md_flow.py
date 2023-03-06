from atomate2.vasp.jobs.core import MDMaker
from jobflow import Flow, Maker
from mpmorph.jobs.core import CHGNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromCalc, PVFromCHGNet, PVFromVasp
from mpmorph.jobs.tasks.chgnet_input import CHGNetMDInputs

EQUILIBRATE_VOLUME_FLOW = "EQUILIBRATE_VOLUME_FLOW"
CHGNET_MD = "CHGNET_MD"

def get_md_flow_chgnet(structure, temp, steps_prod, converge_first = True, initial_vol_scale = 1, use_device='cpu'):
    inputs_prod = CHGNetMDInputs(
        temperature=temp,
        steps=steps_prod,
        use_device=use_device
    )
    inputs_pv = CHGNetMDInputs(
        temperature=temp,
        steps=1000,
        use_device=use_device
    )    
    pv_md_maker = PVFromCHGNet(parameters=inputs_pv)
    chgnet_maker = CHGNetMDMaker(parameters = inputs_prod)
    return _get_md_flow(
        pv_md_maker=pv_md_maker,
        production_md_maker=chgnet_maker,
        structure=structure,
        converge_first=converge_first,
        initial_vol_scale=initial_vol_scale,
    )    


def get_equil_vol_flow_chgnet(structure, temp, steps):
    inputs = CHGNetMDInputs(
        temperature=temp,
        steps=steps
    )

    pv_md_maker = PVFromCHGNet(parameters=inputs)
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)
    equil_vol_job = eq_vol_maker.make(structure)
    flow = Flow([equil_vol_job], output=equil_vol_job.output, name=EQUILIBRATE_VOLUME_FLOW)
    return flow


def _get_md_flow(pv_md_maker, production_md_maker, structure, converge_first, 
initial_vol_scale):

    struct = structure.copy()
    if initial_vol_scale is not None:
        struct.scale_lattice(struct.volume * initial_vol_scale)

    if converge_first:
        return _get_converge_flow(struct, pv_md_maker, production_md_maker)
    else:
        return Flow([production_md_maker.make(struct)], name='flow')

def _get_converge_flow(structure: Structure, pv_md_maker: PVFromCalc, production_run_maker: Maker):
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)

    equil_vol_job = eq_vol_maker.make(structure)

    final_md_job = production_run_maker.make(equil_vol_job.output)

    flow = Flow([equil_vol_job, final_md_job], output=final_md_job.output, name=CHGNET_MD)

    return flow
