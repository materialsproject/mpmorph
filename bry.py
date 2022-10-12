from atomate2.vasp.jobs.core import MDMaker, MDSetGenerator
from jobflow import Flow, job, Response


#pseudo code

#functions and jobs:

#BirchMurnaghan_EOS (directly from MPMorph)
def fit_BirchMurnaghanPV_EOS(p_v):
    # Borrows somewhat from pymatgen/io/abinitio/EOS
    # Initial guesses for the parameters
    from scipy.optimize import leastsq
    eqs = np.polyfit(p_v[:, 1], p_v[:, 0], 2)
    V0 = np.mean(p_v[:, 1])  # still use mean to ensure we are at reasonable volumes
    B0 = -1 * (2 * eqs[0] * V0 ** 2 + eqs[1] * V0)
    B0p = 4.0
    initial_params = (V0, B0, B0p)
    Error = lambda params, x, y: BirchMurnaghanPV_EOS(x, params) - y
    found_params, check = leastsq(Error, initial_params, args=(p_v[:, 1], p_v[:, 0]))
    if check not in [1, 2, 3, 4]:
        raise ValueError("fitting not converged")    
    else:
        return found_params
    

#PV Rescale job
@job
def PV_rescale(md_jobs):
    pv_pairs = np.array([job.output['pressure_volume'] for job in md_jobs]) #(Obtained from MPMorph) it should be tuples
    pv_pairs = np.flip(pv_pairs, axis=1) 
    pv_pairs = np.flip(pv_pairs[pv_pairs[:, 1].argsort()], axis=0)

    params = fit_BirchMurnaghanPV_EOS(pv_pairs)
    equil_volume = params[0]
    
    struc.scale_lattice(equil_volume)
    
    return struc


#flow

#obtains all the structures
rescaler = [0.8, 1, 1.2]
structures = [structure.copy() for factor in rescaler]
for index, factor in enumerate(rescaler):
    structures[index].scale_lattice(structure.volume * factor)
    
#runs an aimd calc for each structure
md_jobs = []
for i, struc in enumerate(structure):
    md_job= MDMaker(name = 'unique identifier_struc_scaler').make(structure)
    md_job.maker.input_set_generator.user_incar_settings["NSW"] = 10  
    md_job.StorePressureVolumeDatainTaskDoc #
    md_jobs.append(md_job)

PV_rescale_job = PV_rescale([job.output for job in md_jobs])
production_run_md_job = MDMaker(name = 'unique identifier_struc_scaler').make(PV_rescale_job.output)

