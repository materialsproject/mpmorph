from distutils.log import ERROR
from doctest import master
from operator import eq
from atomate2.vasp.jobs.core import MDMaker
from jobflow import Flow, job, Response
import numpy as np

from mpmorph.runners import rescale_volume

#flow

def get_converge_flow(structure):
    #obtains all the structures
    scale_factors = [0.8, 1, 1.2]
    md_jobs = [get_md_job(structure,factor) for factor in scale_factors]

    #Dynamic job
    PV_rescale_job = PV_rescale(structure, [job.output for job in md_jobs])
    production_run_md_job = MDMaker().make(PV_rescale_job.output)
    production_run_md_job.maker.input_set_generator.user_incar_settings["NSW"] = NSW  #runs it for short timestep (designed for debugging)

    flow = Flow([*md_jobs, PV_rescale_job, production_run_md_job])
    
    return flow 

#PV Rescale job
@job
def PV_rescale(structure, md_jobs_outputs):
    if len(md_jobs_outputs) > MAX_MD_JOBS:
        raise RuntimeError("YOU'RE DOING SOMETHING WRONG")

    final_structure = structure.copy()

    volume_data = [task_doc_to_volume_data(doc) for doc in md_jobs_outputs]
    pressure_data = [task_doc_to_pressure_data(doc) for doc in md_jobs_outputs]

    max_vol_scale = max(volume_data)
    min_vol_scale = min(volume_data)


    params = rescale_volume.fit_BirchMurnaghanPV_EOS(zip(volume_data,pressure_data))
    equil_volume = params[0]
    
    if equil_volume < max_vol_scale and equil_volume > min_vol_scale:
        return final_structure.scale_lattice(equil_volume)

    elif equil_volume > max_vol_scale: 
        off_vol = get_new_max_volume(equil_volume)

    elif equil_volume < min_vol_scale: 
        off_vol = get_new_min_volume(equil_volume)

    new_job = get_md_job(structure, off_vol)

    md_jobs_outputs.append(new_job.output)
    PV_rescale_job = PV_rescale(structure, md_jobs_outputs)

    flow = Flow([new_job,PV_rescale_job])

    return Response(replace = flow)


def get_new_max_volume(equil_guess, original_structure):
    return equil_guess/original_structure.volume + OFFSET

def get_new_min_volume(equil_guess, original_structure):
    return equil_guess/original_structure.volume - OFFSET

def task_doc_to_volume_data(task_doc): #TODO
    volume = task_doc.calcs_reversed[-1].output.ionic_steps[-1].structure.lattice.volume
    return volume

def task_doc_to_pressure_data(task_doc): #TODO
    stress_tensor = task_doc.calcs_reversed[-1].output.ionic_steps[-1].stress
    pressure = 1/3 * np.trace(stress_tensor)
    return pressure

def get_md_job(structure, factor):
    struct = structure.copy()
    structure.scale_lattice(struct.volume * factor)
    md_job= MDMaker().make(struct)
    md_job.maker.input_set_generator.user_incar_settings["NSW"] = NSW  #runs it for short timestep (designed for debugging)
    md_job.metadata.update({"scale_factor": factor})

    return md_job


#Params list
MAX_MD_JOBS = 5 #if you can't converge with five additional calcs you're doing something wrong...
OFFSET = 0.1 #gives it enough room to slosh back
NSW = 10