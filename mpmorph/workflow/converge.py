from fireworks import Firework, Workflow
from mpmorph.fireworks.core import ConvergeFW
import mpmorph.fireworks.powerups as mp_powerup
from atomate.vasp.fireworks.core import MDFW


def get_converge(structure, temperatures, converge_type="density", unconverged=True, **kwargs):
    kwargs["converge_type"]=converge_type
    fw_list = []

    #TODO: Converge property
    # "density", "total_energy", "kinetic_energy"
    if unconverged:
        fw = MDFW(structure=structure, start_temp=temperature, end_temp=temperature, nsteps=nsteps,
               name=name + "run0", vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, copy_vasp_outputs=False, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)
        mp_powerup.add_converge_task(fw, **kwargs)
        fw_list.append(fw)
        #Add final structure to DB after every stage

    #TODO: Production Run
    #Steps desired, chunks
    while prod_steps < prod_target - prod_length:
        fw = MDFW(structure=structure, start_temp=temperature, end_temp=temperature, nsteps=nsteps,
               name=name + "run0", vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, copy_vasp_outputs=False, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)
        fw = mp_powerup.add_cont_structure(fw, vasp_input_set=vasp_input_set)
        fw_list.append(fw)

    pretty_name=structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name = pretty_name + "_diffusion")