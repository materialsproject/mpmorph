from fireworks import explicit_serialize, Firework, Workflow, FireTaskBase, FWAction
from atomate.vasp.fireworks.core import MDFW
from mpmorph.analysis import md_data
import mpmorph.fireworks.powerups as mp_powerup
import numpy as np

@explicit_serialize
class ConvergeTask(FireTaskBase):
    """

    Ensures a structure is converged before production MD run

    """

    required_params = ["structure", "converge_params", "run_specs", "md_params"]
    optional_params = []

    def run_task(self, fw_spec):
        #Check convergence of all values in converge_params
        converge_params = self["converge_params"]
        data_keys = ['external', 'kinetic energy EKIN', '% ion-electron', 'ETOTAL']
        outcar_data = md_data.get_MD_data("./OUTCAR", search_keys=data_keys)
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
                mp_powerup.add_converge_task(fw, **_spawner_args)
                mp_powerup.add_pass_structure(fw)
                wf = Workflow(fw)
                return FWAction(additions=wf)
        else:
            return FWAction()