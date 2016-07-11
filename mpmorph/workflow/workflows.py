from fireworks import Workflow, Firework
from matmethods.vasp.fireworks.core import MDFW
from mpmorph.workflow.mdtasks import SpawnMDFWTask, CopyCalsHome
from mpmorph.runners.amorphous_maker import AmorphousMaker
from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs
from pymatgen.core.structure import Structure
import os

def get_wf_density(structure, temperature, pressure_threshold=5.0, max_rescales=6, nsteps=2000, wall_time=19200,
                   vasp_input_set=None, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", name="density_finder",
                   optional_MDWF_params=None, amorphous_maker_params=None, copy_calcs=False, calc_home="~/wflows"):

    if copy_calcs:
        if not os.path.exists(calc_home):
            raise ValueError("calc_home must be an existing folder.")
        elif os.path.exists(calc_home+"/"+name):
            raise ValueError("WF name already exists, choose a different name.")
        else:
            calc_home = os.path.join(calc_home,name)
            os.mkdir(calc_home)

    # If structure is in fact just a composition, create a random packed Structure!
    if not isinstance(structure, Structure) and isinstance(structure, dict):
        if not amorphous_maker_params:
            raise ValueError("amorphous_maker_params must be defined!")
        glass = AmorphousMaker(structure, **amorphous_maker_params)
        structure = glass.random_packed_structure

    optional_MDWF_params = optional_MDWF_params or {}
    override_default_vasp_params = {'user_incar_settings' : {"ISIF": 1}}

    fw1 = MDFW(structure=structure, start_temp=temperature, end_temp=temperature, nsteps=nsteps,
               name=name+"run0", vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)
    t = [CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files =["XDATCAR", "OSZICAR"])]

    if copy_calcs:
         t.append(CopyCalsHome(calc_home=calc_home, run_name="run0"))
    t.append(SpawnMDFWTask(pressure_threshold=pressure_threshold, max_rescales=max_rescales,
                       wall_time=wall_time, vasp_cmd=vasp_cmd, db_file=db_file,
                       copy_calcs=copy_calcs, calc_home=calc_home,
                       spawn_count=0))

    fw2 = Firework(t, parents=[fw1], name=name+"_initial_spawn")
    return Workflow([fw1, fw2], name=name+"_WF")