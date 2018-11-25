import os
import numpy as np
from fireworks import explicit_serialize, FireTaskBase, FWAction, Workflow
from pymatgen.io.vasp import Poscar
from pymatgen import Structure, Lattice
from mpmorph.analysis import md_data


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
            for isite, site in enumerate(spec_structure):
                spec_structure[isite] = site.to_unit_cell
            new_lattice = Lattice(_poscar.structure.volume ** (1 / 3.0) * np.eye(3))
            new_fract = new_lattice.get_fractional_coords(spec_structure.cart_coords)
            new_s = Structure(new_lattice, spec_structure.species, new_fract,
                              charge=spec_structure.charge,
                              site_properties=spec_structure.site_properties)
            _structure = new_s.as_dict()

        with open("sst_out", 'w') as f:
            f.write(str(_structure))
            f.close()
        return FWAction(update_spec={"structure": _structure})


@explicit_serialize
class PassPVTask(FireTaskBase):

    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        pressure_volume = fw_spec.get('pressure_volume', [])

        # get volume
        osw = list(os.walk("."))[0]
        files = []
        for file_name in osw[2]:
            if "CONTCAR" in file_name:
                files.append(file_name)
        _poscar = Poscar.from_file(filename=files[-1], check_for_POTCAR=True, read_velocities=True)
        volume = _poscar.structure.volume

        # get pressure
        search_keys = ['external']
        outcar_data = md_data.get_MD_data("./OUTCAR.gz", search_keys=search_keys)

        _data = np.transpose(outcar_data)[0]
        pressure = np.mean(_data[int(0.5 * (len(_data) - 1)):])

        pressure_volume.append((volume, pressure))
        return FWAction(mod_spec={'_push_all': {'pressure_volume': pressure_volume}})
