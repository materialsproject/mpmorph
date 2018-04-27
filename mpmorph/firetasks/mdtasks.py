from fireworks import explicit_serialize, Workflow, FireTaskBase, FWAction
from mpmorph.runners.rescale_volume import RescaleVolume
from mpmorph.util import recursive_update
from mpmorph.analysis import md_data
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import numpy as np
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor

@explicit_serialize
class DLSVPRescaling(FireTaskBase):
    """
    After First MD run, manipulate the final structure according to:
    1) Rescale according to DSLVolumePredictor
    2) Run an optimize FW.
    3) Pass Structure
    """
    required_params = []
    optional_params = []

    def run_task(self, fw_spec):
        structure = Structure.from_dict(fw_spec["structure"])
        dlsvp = DLSVolumePredictor()
        dlsvp_structure = dlsvp.predict(structure)
        return FWAction(update_spec={"structure": dlsvp_structure})

@explicit_serialize
class ConvergeTask(FireTaskBase):
    """

    Ensures a structure is converged before production MD run

    """

    required_params = ["converge_params", "run_specs", "md_params"]
    optional_params = ["rescale_params"]

    def run_task(self, fw_spec):
        # Load Structure from Poscar
        _poscar = Poscar.from_file("CONTCAR.gz")
        structure = _poscar.structure

        #Get convergence parameters from spec
        converge_params = self["converge_params"]
        avg_fraction = converge_params.get("avg_fraction", 0.5)
        convergence_vars = dict(converge_params["converge_type"])
        if "total energy" not in convergence_vars.keys():
            convergence_vars["total_energy"]=0.0005
        rescale_params = self.get("rescale_params", {})

        #Load Data from OUTCAR
        key_map = {'density': 'external', 'kinetic energy': 'kinetic energy EKIN', 'ionic': '% ion-electron', 'total energy': 'ETOTAL'}
        outcar_data = md_data.get_MD_data("./OUTCAR.gz", search_keys=list(key_map.values()))


        # Check for convergence
        converged={}
        _index = list(key_map.values()).index(key_map["density"])
        _data = np.transpose(outcar_data)[_index].copy()
        pressure = np.mean(_data[int(avg_fraction * (len(_data) - 1)):])
        if "density" in convergence_vars.keys():
            if np.abs(pressure) >= convergence_vars["density"]:
                converged["density"] = False
            else:
                converged["density"] = True


        if "kinetic energy" in convergence_vars.keys():
            _index = list(key_map.values()).index(key_map["kinetic_energy"])
            energy = np.transpose(outcar_data)[_index].copy()
            norm_energy = (energy / structure.num_sites) / np.mean(energy / structure.num_sites) - 1
            if np.abs(np.mean(norm_energy[-500:])-np.mean(norm_energy)) > convergence_vars["kinetic energy"]:
                converged["kinetic_energy"] = False
            else:
                converged["kinetic_energy"] = True

        _index = list(key_map.values()).index(key_map["total_energy"])
        energy = np.transpose(outcar_data)[_index].copy()
        norm_energy = (energy / structure.num_sites) / np.mean(energy / structure.num_sites) - 1
        if np.abs(np.mean(norm_energy[-500:]) - np.mean(norm_energy)) > convergence_vars["total_energy"]:
            converged["total_energy"] = False
        else:
            converged["total_energy"] = True

        # Spawn Additional Fireworks
        if not all([item[1] for item in converged.items()]):
            spawn_count = converge_params["spawn_count"]
            max_spawns = converge_params["max_rescales"]
            if spawn_count >= max_spawns:
                return FWAction(defuse_children=True)
            else:
                run_specs = self["run_specs"]
                md_params = self["md_params"]
                optional_params = self["optional_fw_params"]
                if not converged["density"]:
                    rescale_args = {"initial_pressure": pressure * 1000, "initial_temperature": 1, "beta": 0.0000005}
                    rescale_args = recursive_update(rescale_args, rescale_params)

                    # Spawn fw
                    fw = MDFW(structure, name="density_run" + str(spawn_count + 1), previous_structure=False,
                              insert_db=False, **run_specs, **md_params, **optional_params)
                    converge_params["spawn_count"] += 1
                    _spawner_args = {"converge_params": converge_params, "rescale_params": rescale_params,
                                     "run_specs": run_specs, "md_params": md_params,
                                     "optional_fw_params": optional_params}
                    fw = powerups.add_rescale_volume(fw, **rescale_args)
                    fw = powerups.add_converge_task(fw, **_spawner_args)
                else:
                    fw = MDFW(structure, name="energy_run" + str(spawn_count + 1), previous_structure=True,
                              insert_db=False, **run_specs, **md_params, **optional_params)
                    converge_params["spawn_count"] += 1
                    _spawner_args = {"converge_params": converge_params, "rescale_params": rescale_params,
                                     "run_specs": run_specs, "md_params": md_params,
                                     "optional_fw_params": optional_params}
                    fw = powerups.add_converge_task(fw, **_spawner_args)
                wf = Workflow([fw])
                return FWAction(detours=wf, stored_data={'pressure': pressure})
        else:
            return FWAction(stored_data={'pressure': pressure, 'density_calculated': True})


@explicit_serialize
class RescaleVolumeTask(FireTaskBase):
    """
    Volume rescaling
    """
    required_params = ["initial_temperature", "initial_pressure"]
    optional_params = ["target_pressure", "target_temperature", "target_pressure", "alpha", "beta"]

    def run_task(self, fw_spec):
        # Initialize volume correction object with last structure from last_run
        initial_temperature = self["initial_temperature"]
        initial_pressure = self["initial_pressure"]
        target_temperature = self.get("target_temperature", initial_temperature)
        target_pressure = self.get("target_pressure", 0.0)
        alpha = self.get("alpha", 10e-6)
        beta = self.get("beta", 10e-7)
        corr_vol = RescaleVolume.of_poscar(poscar_path="./POSCAR", initial_temperature=initial_temperature,
                                           initial_pressure=initial_pressure,
                                           target_pressure=target_pressure,
                                           target_temperature=target_temperature, alpha=alpha, beta=beta)
        # Rescale volume based on temperature difference first. Const T will return no volume change:
        corr_vol.by_thermo(scale='temperature')
        # TO DB ("Rescaled volume due to delta T: ", corr_vol.structure.volume)
        # Rescale volume based on pressure difference:
        corr_vol.by_thermo(scale='pressure')
        # TO DB ("Rescaled volume due to delta P: ", corr_vol.structure.volume)
        corr_vol.poscar.write_file("./POSCAR")
        # Pass the rescaled volume to Poscar
        return FWAction(stored_data=corr_vol.structure.as_dict())

from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW
