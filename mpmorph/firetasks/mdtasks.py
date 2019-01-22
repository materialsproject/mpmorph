from scipy import stats
import numpy as np
from fireworks import explicit_serialize, Workflow, FireTaskBase, FWAction
from pymatgen import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor
from mpmorph.runners.rescale_volume import RescaleVolume, fit_BirchMurnaghanPV_EOS
from mpmorph.util import recursive_update
from mpmorph.analysis import md_data


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
        with open("sst_out", 'w') as f:
            f.write(str(fw_spec["structure"]))
            f.close()
        structure = fw_spec["structure"]
        dlsvp = DLSVolumePredictor()
        dlsvp_volume = dlsvp.predict(structure)
        structure.scale_lattice(dlsvp_volume)
        return FWAction(update_spec={"structure": structure})


@explicit_serialize
class DiffusionTask(FireTaskBase):
    required_params = ['temperatures', 'max_steps', 'target_steps',
                       'trajectory_to_db', 'notes']
    optional_params = []

    def run_task(self, fw_spec):
        from mpmorph.workflow.converge import get_converge_new

        s = Structure.from_file('CONTCAR.gz')
        fws = []
        for t in self['temperatures']:
            fws.extend(get_converge_new(s, int(t), max_steps=self['max_steps'],
                                        target_steps=self['target_steps'],
                                        trajectory_to_db=self['trajectory_to_db'],
                                        notes=self['notes']))
        wf = Workflow(fws)
        return FWAction(detours=wf)


@explicit_serialize
class ConvergeTask(FireTaskBase):
    """

    Ensures a structure is converged before production MD run

    """

    required_params = ["converge_params", "run_specs", "md_params"]
    optional_params = ["rescale_params"]

    def run_task(self, fw_spec):
        from mpmorph.fireworks import powerups
        from mpmorph.fireworks.core import MDFW

        # Load Structure from Poscar
        _poscar = Poscar.from_file("CONTCAR.gz")
        structure = _poscar.structure

        # Get convergence parameters from spec
        converge_params = self["converge_params"]
        avg_fraction = converge_params.get("avg_fraction", 0.5)
        convergence_vars = dict(converge_params["converge_type"])
        if "ionic" not in convergence_vars.keys():
            convergence_vars["ionic"]=0.0005
        rescale_params = self.get("rescale_params", {})

        # Load Data from OUTCAR
        search_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL']
        key_map = {'density': 'external', 'kinetic energy': 'kinetic energy EKIN', 'ionic': '% ion-electron',
                   'total energy': 'ETOTAL'}
        outcar_data = md_data.get_MD_data("./OUTCAR.gz", search_keys=search_keys)

        # Check for convergence
        converged={}
        _index = search_keys.index(key_map["density"])
        _data = np.transpose(outcar_data)[_index].copy()
        pressure = np.mean(_data[int(avg_fraction * (len(_data) - 1)):])
        if "density" in convergence_vars.keys():
            if np.abs(pressure) >= convergence_vars["density"]:
                converged["density"] = False
            else:
                converged["density"] = True

        if "kinetic energy" in convergence_vars.keys():
            _index = search_keys.index(key_map["kinetic energy"])
            energy = np.transpose(outcar_data)[_index].copy()
            norm_energy = (energy / structure.num_sites) / np.mean(energy / structure.num_sites) - 1
            if np.abs(np.mean(norm_energy[-500:])-np.mean(norm_energy)) > convergence_vars["kinetic energy"]:
                converged["kinetic energy"] = False
            else:
                converged["kinetic energy"] = True

        _index = search_keys.index(key_map["ionic"])
        energy = np.transpose(outcar_data)[_index].copy()
        norm_energy = (energy / structure.num_sites) / np.mean(energy / structure.num_sites) - 1
        norm_energy_change = np.abs(np.mean(norm_energy[-500:]) - np.mean(norm_energy))
        if norm_energy_change > convergence_vars["ionic"]:
            converged["ionic"] = False
        else:
            converged["ionic"] = True

        # Spawn Additional Fireworks
        if not all([item[1] for item in converged.items()]):
            spawn_count = converge_params["spawn_count"]
            max_spawns = converge_params["max_rescales"]

            run_specs = self["run_specs"]
            md_params = self["md_params"]
            optional_params = self["optional_fw_params"]
            if spawn_count >= max_spawns:
                return FWAction(defuse_children=True)
            elif not converged.get("density", True):
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
                fw = MDFW(structure, name="energy_run" + str(spawn_count + 1), previous_structure=False,
                          insert_db=False, **run_specs, **md_params, **optional_params)
                converge_params["spawn_count"] += 1
                _spawner_args = {"converge_params": converge_params, "rescale_params": rescale_params,
                                 "run_specs": run_specs, "md_params": md_params,
                                 "optional_fw_params": optional_params}
                fw = powerups.add_converge_task(fw, **_spawner_args)
            wf = Workflow([fw])
            return FWAction(detours=wf, stored_data={'pressure': pressure,
                                                     'energy': norm_energy_change})
        else:
            return FWAction(stored_data={'pressure': pressure,
                                         'energy': norm_energy_change,
                                         'density_calculated': True})


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


@explicit_serialize
class PVRescaleTask(FireTaskBase):
    """
    Rescale based on fitting pressure vs volume to Birch-Murnaghan EOS
    """

    required_params = []
    optional_params = ['rescale_type']

    def run_task(self, fw_spec):
        rescale_type = self.get('rescale_type', 'BirchMurnaghan_EOS')

        if rescale_type == 'BirchMurnaghan_EOS':
            pv_pairs = np.array(fw_spec["pressure_volume"])
            pv_pairs = np.flip(pv_pairs, axis=1)
            pv_pairs = np.flip(pv_pairs[pv_pairs[:, 1].argsort()], axis=0)

            params = fit_BirchMurnaghanPV_EOS(pv_pairs)
            equil_volume = params[0]
        else:
            pvs = fw_spec["pressure_volume"]

            p = [item[1] for item in pvs]
            v = [item[0] for item in pvs]

            slope, intercept, r_value, p_value, std_err = stats.linregress(v, p)
            equil_volume = -intercept/slope

        poscar = Poscar.from_file("./POSCAR")
        poscar.structure.scale_lattice(equil_volume)
        poscar.write_file("./POSCAR")

        return FWAction()
