from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from pymatgen import Structure

@explicit_serialize
class PreviousStructureTask(FireTaskBase):

    def run_task(self, fw_spec):
        #get last structure from fw_spec
        poscar_dicts = fw_spec['poscar']
        _poscar = Poscar.from_dict(poscar_dicts[-1])
        _poscar.write_file("POSCAR")
        return FWAction()
