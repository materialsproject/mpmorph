from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from mpmorph.runners.rescale_volume import RescaleVolume
from pymatgen.io.vasp import Poscar
from mpmorph.analysis import md_data
import numpy as np
from mpmorph.util import recursive_update

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

        #Check convergence of all values in converge_params
        converge_params = self["converge_params"]
        rescale_params = self.get("rescale_params", {})
        data_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL']
        key_map = {'density': 'external', 'kinetic energy': 'kinetic energy EKIN', 'ionic': '% ion-electron', 'total energy': 'ETOTAL'}
        outcar_data = md_data.get_MD_data("./OUTCAR.gz", search_keys=data_keys)
        convergence_vars = converge_params["converge_type"]
        converged = False
        pressure = 0

        for converge_tuple in convergence_vars:
            converge_key, threshold = converge_tuple

            avg_fraction = converge_params.get("avg_fraction", 0.5)

            _index = data_keys.index(key_map[converge_key])
            _data = np.transpose(outcar_data)[_index]
            avg_val = np.mean(_data[int(avg_fraction*(len(_data)-1)):])
            pressure = avg_val

            # Pressure Convergence
            if key_map[converge_key] == 'external':
                converged = True if np.abs(avg_val) <= threshold else False

            # Energy Convergence
            if key_map[converge_key] in ['kinetic energy EKIN', 'ETOTAL']:
                continue

        if not converged:
            spawn_count = converge_params["spawn_count"]
            max_spawns = converge_params["max_rescales"]
            if spawn_count >= max_spawns:
                return FWAction(defuse_children=True)
            else:
                run_specs = self["run_specs"]
                md_params = self["md_params"]
                optional_params = self["optional_fw_params"]

                rescale_args = {"initial_pressure": pressure*1000, "initial_temperature": 1, "beta": 0.000002}
                rescale_args = recursive_update(rescale_args, rescale_params)

                #Spawn fw
                fw = MDFW(structure, name="run"+ str(spawn_count+1), previous_structure=False, insert_db=False, **run_specs, **md_params, **optional_params)
                converge_params["spawn_count"] += 1
                _spawner_args = {"converge_params": converge_params, "rescale_params": rescale_params, "run_specs": run_specs, "md_params":md_params, "optional_fw_params": optional_params}
                fw = powerups.add_rescale_volume(fw, **rescale_args)
                fw = powerups.add_converge_task(fw, **_spawner_args)
                wf = Workflow([fw])
                return FWAction(detours=wf)
        else:
            return FWAction()

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
