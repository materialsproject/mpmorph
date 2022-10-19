from atomate2.vasp.jobs.core import MDMaker
from jobflow import Flow, Maker
from mpmorph.flows.scale_volume import md_to_volume_flow
from mpmorph.jobs.core import M3GNetMDMaker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure
from mpmorph.jobs.extract_pv_m3gnet import ExtractPVDataM3gnet

from mpmorph.jobs.extract_pv_vasp import ExtractPVDataFromVASPMDMaker

def get_converge_flow_vasp(structure, temperature):
    md_maker = MDMaker()
    NSW = 10
    md_maker.input_set_generator.user_incar_settings["NSW"] = NSW  #runs it for short timestep (designed for debugging)
    pv_extract_maker = ExtractPVDataFromVASPMDMaker()
    return get_converge_flow(structure, md_maker, pv_extract_maker)

def get_converge_flow_m3gnet(structure, temp):
    md_maker = M3GNetMDMaker(
        temperature=temp
    )
    pv_extract_maker = ExtractPVDataM3gnet()
    return get_converge_flow(structure, md_maker, pv_extract_maker)

def get_converge_flow(structure: Structure, md_maker: Maker, pv_extract_maker: Maker):
    # TODO: Make this initial range tunable - in some cases, these may lead
    # to unwanted(?) phase changes
    initial_scale_factors = [0.8, 1, 1.2]
    
    initial_vol_search_jobs = [
        md_to_volume_flow(
            structure,
            factor,
            md_maker,
            pv_extract_maker
        ) for factor in initial_scale_factors
    ]

    eq_vol_maker = EquilibriumVolumeSearchMaker(
        md_maker = md_maker,
        pv_maker = pv_extract_maker
    )

    equil_vol_job = eq_vol_maker.make(structure, [job.output for job in initial_vol_search_jobs])
    
    final_md_job = md_maker.make(equil_vol_job.output)

    flow = Flow([*initial_vol_search_jobs, equil_vol_job, final_md_job])
    
    return flow 
