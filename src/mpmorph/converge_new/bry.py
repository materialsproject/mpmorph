from distutils.log import ERROR
from operator import eq
from atomate2.vasp.jobs.core import MDMaker, MDSetGenerator
from jobflow import Flow, job, Response

from mpmorph.runners import rescale_volume

#flow

def get_converge_flow(structure):
    #obtains all the structures
    scale_factors = [0.8, 1, 1.2]
    md_jobs = [get_md_job(structure,factor) for factor in scale_factors]

    #Dynamic job
    PV_rescale_job = PV_rescale(structure, [job.outputs for job in md_jobs])
    production_run_md_job = MDMaker(name = 'unique identifier_struc_scaler').make(PV_rescale_job.output)
    
    flow = Flow([md_jobs, PV_rescale_job, production_run_md_job])
    
    return flow 

#PV Rescale job
@job
def PV_rescale(structure, md_jobs_outputs):
    if len(md_jobs_outputs) > MAX_MD_JOBS:
        raise RuntimeError("YOU'RE DOING SOMETHING WRONG")

    final_structure = structure.copy()

    pv_pairs = np.array([job['pressure_volume'] for job in md_jobs_outputs]) #(Obtained from MPMorph) it should be tuples
    pv_pairs = np.flip(pv_pairs, axis=1) 
    pv_pairs = np.flip(pv_pairs[pv_pairs[:, 1].argsort()], axis=0)

    max_vol_scale = 
    min_vol_scale =


    params = fit_BirchMurnaghanPV_EOS(pv_pairs)
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

MAX_MD_JOBS = 5 #if you can't converge with five additional calcs you're doing something wrong...
OFFSET = 0.1 #gives it enough room to slosh back
def get_new_max_volume(equil_guess, original_structure):
    return equil_guess/original_structure.volume + OFFSET

def get_new_min_volume(equil_guess, original_structure):
    return equil_guess/original_structure.volume - OFFSET

def task_doc_to_volume_data(task_doc): #TODO
    volume = task_doc['output']['output']['stress']

def task_doc_to_pressure_data(task_doc): #TODO
    stress_tensor = task_doc['output']['output']['stress']
    pressure = trace(stress_tensor)
    return pressure

def get_md_job(structure, factor):
    struct = structure.copy()
    structure.scale_lattice(struct.volume * factor)
    md_job= MDMaker().make(struct)
    md_job.maker.input_set_generator.user_incar_settings["NSW"] = 10  
    md_job.metadata.update({"scale_factor": factor})
    md_job.StorePressureVolumeDatainTaskDoc #

    return md_job

