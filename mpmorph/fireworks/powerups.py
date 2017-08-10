from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassResult
from pymatgen.io.vasp import Poscar
from mpmorph.firetasks.mdtasks import ConvergeTask
from mpmorph.firetasks.util import PreviousStructureTask

def add_converge_task(fw, **kwargs):
    #Load Structure from Poscar
    _poscar = Poscar.from_file("CONTCAR")
    structure = _poscar.structure

    spawner_task = ConvergeTask(structure, **kwargs)
    fw.tasks.append(spawner_task)
    return fw

def add_cont_structure(fw, position):
    t = PreviousStructureTask()
    fw.tasks.append(position, t)
    return fw

def add_pass_structure(fw, **kwargs):
    pass_structure = PassResult(pass_dict={"structure": ">>structures.-1"}, parse_class="pymatgen.io.vasp.Vasprun",
                                parse_kwargs={"filename": "vasprun.xml", "parse_dos": False, "parse_eigen": False})
    fw.tasks.append(pass_structure)
    return fw