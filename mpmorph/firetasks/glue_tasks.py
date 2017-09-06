from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from pymatgen import Structure

@explicit_serialize
class PreviousStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        #get last structure from fw_spec
        structure_dict = fw_spec["structure"]
        _structure = Structure.from_dict(structure_dict)
        _poscar = Poscar(_structure)
        _poscar.write_file("POSCAR")
        return FWAction()

@explicit_serialize
class SaveStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        _poscar = Poscar.from_file(filename="CONTCAR.gz", check_for_POTCAR=True, read_velocities=True)
        _structure = _poscar.structure.as_dict()

        return FWAction(mod_spec=[{ "_push" : {"structure": _structure}}])