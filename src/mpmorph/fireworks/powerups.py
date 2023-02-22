from mpmorph.firetasks.dbtasks import TrajectoryDBTask, VaspMDToDb
from mpmorph.firetasks.glue_tasks import (
    PassPVTask,
    PreviousStructureTask,
    SaveStructureTask,
)
from mpmorph.firetasks.mdtasks import (
    ConvergeTask,
    DiffusionTask,
    PVRescaleTask,
    RescaleVolumeTask,
)


def add_diffusion_task(fw, **kwargs):
    spawner_task = DiffusionTask(**kwargs)
    fw.tasks.append(spawner_task)
    return fw


def add_converge_task(fw, **kwargs):
    """
    This powerup adds the convergence task onto a MD firework, which turns a workflow into a dynamic workflow.
    This firetask will check to ensure the specified parameters (Usually pressure and Energy) are converged within the
    specified thresholds.
    :param fw:
    :param kwargs:
    :return:
    """
    spawner_task = ConvergeTask(**kwargs)
    fw.tasks.append(spawner_task)
    return fw


def aggregate_trajectory(fw, **kwargs):
    """
    This firetask will add a task which converts a series of MD runs into a trajectory object
    :param fw:
    :param kwargs:
    :return:
    """
    fw.tasks.append(TrajectoryDBTask(**kwargs))
    return fw


def add_cont_structure(fw):
    prev_struct_task = PreviousStructureTask()
    insert_i = 2
    for i, task in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i
            break
    fw.tasks.insert(insert_i, prev_struct_task)
    return fw


def add_pass_structure(fw, **kwargs):
    save_struct_task = SaveStructureTask(**kwargs)
    fw.tasks.append(save_struct_task)
    return fw


def add_pass_pv(fw, **kwargs):
    pass_pv_task = PassPVTask(**kwargs)
    fw.tasks.append(pass_pv_task)
    return fw


def add_pv_volume_rescale(fw):
    insert_i = 2
    for i, task in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i
            break

    fw.tasks.insert(insert_i, PVRescaleTask())
    return fw


def add_rescale_volume(fw, **kwargs):
    rsv_task = RescaleVolumeTask(**kwargs)
    insert_i = 2
    for i, task in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i
            break

    fw.tasks.insert(insert_i, rsv_task)
    return fw


def replace_pass_structure(fw, **kwargs):
    # look for rescale_volume task
    replaced = False
    fw_dict = fw.to_dict()
    for i in range(len(fw_dict["spec"]["_tasks"])):
        if (
            fw_dict["spec"]["_tasks"][i]["_fw_name"]
            == "{{mpmorph.firetasks.glue_tasks.SaveStructureTask}}"
        ):
            del fw_dict["spec"]["_tasks"][i]["_fw_name"]
            fw.tasks[i] = SaveStructureTask(**kwargs)
            replaced = True
            break
    # TODO: Replace with real error handling
    if replaced == False:
        print("error, no SaveStructureTask to replace")
        return

    return fw


def replace_vaspmdtodb(fw):
    # look for vaspdb task
    replaced = False
    fw_dict = fw.to_dict()
    for i in range(len(fw_dict["spec"]["_tasks"])):
        if (
            fw_dict["spec"]["_tasks"][i]["_fw_name"]
            == "{{atomate.vasp.firetasks.parse_outputs.VaspToDb}}"
        ):
            del fw_dict["spec"]["_tasks"][i]["_fw_name"]
            fw.tasks[i] = VaspMDToDb(**fw_dict["spec"]["_tasks"][i])
            replaced = True
            break
    # TODO: Replace with real error handling
    if replaced == False:
        print("error, no vasptodb to replace")
        return

    return fw
