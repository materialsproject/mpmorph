from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from pymatgen.io.vasp import Poscar
from atomate.vasp.fireworks.core import MDFW
from mpmorph.analysis import md_data
import numpy as np

@explicit_serialize
class ConvergeTask(FireTaskBase):
    """

    Ensures a structure is converged before production MD run

    """

    required_params = ["converge_params", "run_specs", "md_params"]
    optional_params = []

    def run_task(self, fw_spec):
        # Load Structure from Poscar

        _poscar = Poscar.from_file("CONTCAR.gz")
        structure = _poscar.structure

        #Check convergence of all values in converge_params
        converge_params = self["converge_params"]
        data_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL']
        outcar_data = md_data.get_MD_data("./OUTCAR.gz", search_keys=data_keys)
        convergence_vars = converge_params["converge_type"]
        converged = False

        for converge_tuple in convergence_vars:
            converge_key, threshold = converge_tuple
            avg_fraction = converge_params["avg_fraction"]

            _index = data_keys.index(converge_key)
            _data = np.transpose(outcar_data)[_index]
            avg_val = np.mean(_data[int(avg_fraction*(len(_data)-1)):])

            # Pressure Convergence
            if converge_key == 'external':
                if avg_val <= threshold:
                    converged = True
                else:
                    converged = False

            # Energy Convergence
            if converge_key in ['kinetic energy EKIN', 'ETOTAL']:
                continue


        if not converged:
            spawn_count = converge_params["spawn_count"]
            max_spawns = converge_params["max_rescales"]
            if spawn_count >= max_spawns:
                return FWAction(defuse_children=True)
            else:
                structure = self["structure"]
                run_specs = self["run_specs"]
                md_params = self["md_params"]

                #Spawn fw
                fw = MDFW(structure, **run_specs, **md_params)
                _spawner_args = {"converge_params":converge_params, "run_specs":run_specs, "md_params":md_params}
                fw = powerups.replace_vaspmdtodb(fw)
                fw = powerups.add_cont_structure(fw, position=1) #Add after MDFW WriteInputSet to override structure
                fw = powerups.add_converge_task(fw, **_spawner_args)
                fw = powerups.add_pass_structure(fw)
                wf = Workflow(fw)
                return FWAction(additions=wf)
        else:
            return FWAction()

from mpmorph.fireworks import powerups