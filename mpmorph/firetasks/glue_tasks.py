from fireworks import explicit_serialize, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import os


@explicit_serialize
class PreviousStructureTask(FireTaskBase):

    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        structure_dict = fw_spec["structure"]
        _poscar = Poscar(structure_dict)
        _poscar.write_file("POSCAR")
        return FWAction()


@explicit_serialize
class SaveStructureTask(FireTaskBase):

    required_params = []
    optional_params = ["rescale_volume"]

    def run_task(self, fw_spec):
        osw = list(os.walk("."))[0]
        files = []
        for file_name in osw[2]:
            if "CONTCAR" in file_name:
                files.append(file_name)
        _poscar = Poscar.from_file(filename=files[-1], check_for_POTCAR=True, read_velocities=True)
        _structure = _poscar.structure.as_dict()

        if self.get("rescale_volume", False):
            spec_structure = Structure.from_dict(fw_spec["structure"])
            scaling_volume = _poscar.structure.volume()
            spec_structure.scale_lattice(scaling_volume)
            _structure = spec_structure.as_dict()

        with open("sst_out", 'w') as f:
            f.write(str(_structure))
            f.close()
        return FWAction(update_spec={"structure": _structure})
