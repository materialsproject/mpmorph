from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import os

@explicit_serialize
class PreviousStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        #get last structure from fw_spec
        structure_dict = fw_spec["structure"]
        print(len(structure_dict))
        print(type(structure_dict))
        print(type(structure_dict[-1]))
        print(structure_dict)
        _poscar = Poscar(structure_dict)
        _poscar.write_file("POSCAR")
        return FWAction()

@explicit_serialize
class SaveStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        osw = list(os.walk("."))[0]
        files = []
        for file_name in osw[2]:
            if "CONTCAR" in file_name:
                files.append(file_name)
        _poscar = Poscar.from_file(filename=files[-1], check_for_POTCAR=True, read_velocities=True)
        _structure = _poscar.structure.as_dict()
        with open("sst_out", 'w') as f:
            f.write(str(_structure))
            f.close()
        return FWAction(update_spec={"structure": _structure})
        #return FWAction(mod_spec=[{ "_push" : {"structure": _structure}}])
