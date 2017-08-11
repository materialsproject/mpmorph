from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassResult
from fireworks import Firework

from mpmorph.firetasks.mdtasks import ConvergeTask
from mpmorph.firetasks.util import PreviousStructureTask
from mpmorph.firetasks.dbtasks import VaspMDToDb

def add_converge_task(fw, **kwargs):
    spawner_task = ConvergeTask(**kwargs)
    fw.tasks.append(spawner_task)
    return fw

def add_cont_structure(fw, position):
    prev_struct_task = PreviousStructureTask()
    fw.tasks.insert(position, prev_struct_task)
    return fw

def add_pass_structure(fw, **kwargs):
    pass_structure = PassResult(pass_dict={"structure": ">>structures.-1"}, parse_class="pymatgen.io.vasp.Vasprun",
                                parse_kwargs={"filename": "vasprun.xml", "parse_dos": False, "parse_eigen": False})
    fw.tasks.append(pass_structure)
    return fw

def replace_vaspmdtodb(fw, **kwargs):
    #look for vaspdb task
    replaced = False
    fw_dict = fw.to_dict()
    for i in range(len(fw_dict['spec']['_tasks'])):
        print(fw_dict['spec']['_tasks'][i]["_fw_name"])
        if fw_dict['spec']['_tasks'][i]["_fw_name"] == '{{atomate.vasp.firetasks.parse_outputs.VaspToDb}}':
            del fw_dict['spec']['_tasks'][i]["_fw_name"]
            fw.tasks[i] = VaspMDToDb(**fw_dict['spec']['_tasks'][i])
            print(fw.tasks)
            replaced = True
            break
    #TODO: Replace with real error handling
    if replaced == False:
        print("error, no vasptodb to replace")
        return

    return fw
