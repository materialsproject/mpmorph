from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassResult
from mpmorph.firetasks.mdtasks import RescaleVolumeTask
from mpmorph.firetasks.vibtasks import AdsorbateGeneratorTask
from pymatgen.io.vasp.inputs import Poscar
from fireworks import Firework

from mpmorph.firetasks.mdtasks import ConvergeTask
from mpmorph.firetasks.glue_tasks import PreviousStructureTask, SaveStructureTask
from mpmorph.firetasks.dbtasks import VaspMDToDb

def add_converge_task(fw, **kwargs):
    spawner_task = ConvergeTask(**kwargs)
    insert_i = -2
    for (i, task) in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i+1
            break

    #fw.tasks.insert(insert_i, spawner_task)
    fw.tasks.append(spawner_task)
    return fw

def add_cont_structure(fw):
    prev_struct_task = PreviousStructureTask()
    insert_i = 2
    for (i, task) in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i
            break
    fw.tasks.insert(insert_i, prev_struct_task)
    return fw

def add_pass_structure(fw, velocity=True, **kwargs):
    save_struct_task = SaveStructureTask()
    fw.tasks.append(save_struct_task)
    return fw


def add_rescale_volume(fw, **kwargs):
    rsv_task = RescaleVolumeTask(**kwargs)
    insert_i = 2
    for (i, task) in enumerate(fw.tasks):
        if task.fw_name == "{{atomate.vasp.firetasks.run_calc.RunVaspCustodian}}":
            insert_i = i
            break

    fw.tasks.insert(insert_i, rsv_task)
    return fw

def replace_vaspmdtodb(fw, **kwargs):
    #look for vaspdb task
    replaced = False
    fw_dict = fw.to_dict()
    for i in range(len(fw_dict['spec']['_tasks'])):
        #print(fw_dict['spec']['_tasks'][i]["_fw_name"])
        if fw_dict['spec']['_tasks'][i]["_fw_name"] == '{{atomate.vasp.firetasks.parse_outputs.VaspToDb}}':
            del fw_dict['spec']['_tasks'][i]["_fw_name"]
            fw.tasks[i] = VaspMDToDb(**fw_dict['spec']['_tasks'][i])
            #print(fw.tasks)
            replaced = True
            break
    #TODO: Replace with real error handling
    if replaced == False:
        print("error, no vasptodb to replace")
        return

    return fw

def add_adsorbate_task(fw, **kwargs):
    asd_task = AdsorbateGeneratorTask(**kwargs)
    fw.tasks.append(asd_task)
    return fw