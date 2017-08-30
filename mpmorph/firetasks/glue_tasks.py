from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from pymatgen import Structure

@explicit_serialize
class PreviousStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        #get last structure from fw_spec
        poscar_dict = fw_spec['CONTCAR']
        _poscar = Poscar.from_dict(poscar_dict)
        _poscar.write_file("POSCAR")
        return FWAction()

@explicit_serialize
class SaveStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        _poscar = Poscar.from_file(filename="CONTCAR.gz", check_for_POTCAR=True, read_velocities=True)
        poscar_dict = _poscar.as_dict()

        with open("testwrite", 'w') as file:
            file.write(_poscar.get_string())
            file.close()
        return FWAction(mod_spec=[{ "_push" : {"CONTCAR": poscar_dict}}])