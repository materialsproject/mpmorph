from fireworks import Workflow, Firework
from matmethods.vasp.fireworks.core import MDFW
from mpmorph.workflow.mdtasks import GetPressureTask, SpawnMDFWTask
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs

def get_wf_density(structure, temperature, pressure_threshold=5.0, max_rescales=6, nsteps=2000, wall_time=19200,
                       vasp_input_set=None, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", optional_MDWF_params=None):

    optional_MDWF_params = optional_MDWF_params or {}
    name="md density finder"
    override_default_vasp_params = {'user_incar_settings' : {"ISIF": 1}}

    fw1 = MDFW(structure=structure, start_temp=temperature, end_temp=temperature, nsteps=nsteps,
               name=name, vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)

    fw2 = Firework([CopyVaspOutputs(calc_loc=True), GetPressureTask(outcar_path='./'), PassCalcLocs(name=name)] , parents=[fw1])
    fw3 = Firework([SpawnMDFWTask(pressure_threshold=pressure_threshold,
                                  max_rescales=max_rescales,
                                  wall_time=wall_time, vasp_cmd=vasp_cmd, db_file=db_file)], parents=[fw2])
    return Workflow([fw1, fw2, fw3])