from fireworks import Workflow
from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW, OptimizeFW
from mpmorph.util import recursive_update
from copy import deepcopy
import uuid


def get_converge(structure, priority=None, preconverged=False, max_steps=5000, target_steps=40000,
                 spawner_args={}, converge_args={}, prod_args={}, converge_type=("density", 5),
                 vasp_opt=True, parents=None, descriptor='', **kwargs):
    fw_list = []
    # Initial Run and convergence of structure

    run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<",
                              "db_file": ">>db_file<<", "wall_time": 86400},
                "optional_fw_params": {"override_default_vasp_params": {},
                                       "copy_vasp_outputs": False, "spec": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'ISIF': 1, 'LWAVE': False,
                                 'PREC': 'Normal', 'EDIFF': 1E-4}})
    run_args["optional_fw_params"]["spec"]["_queueadapter"] = {"walltime": run_args["run_specs"]["wall_time"]}
    run_args = recursive_update(run_args, converge_args)
    run_args["optional_fw_params"]["spec"]["_priority"] = priority

    _spawner_args = {"converge_params": {"converge_type": [converge_type], "max_rescales": 15, "spawn_count": 1},
                     "rescale_params": {"beta": 5e-6},
                     "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
                     "optional_fw_params": run_args["optional_fw_params"]}
    _spawner_args["md_params"].update({"start_temp": run_args["md_params"]["end_temp"]})
    _spawner_args = recursive_update(_spawner_args, spawner_args)

    if not preconverged:
        if vasp_opt:
            """
            fw1 for AIMD atomic position optimization.
            fw2 for DFT volume optimization.
            fw3 rescale volume for AIMD atomic position optimization until
            pressure is 0.
            """
            fw1 = MDFW(structure=structure, name=descriptor + 'run0', insert_db=False,
                       previous_structure=False, parents=parents, **run_args["md_params"],
                       **run_args["run_specs"], **run_args["optional_fw_params"])

            # OptimizeFW does not take wall_time
            optimize_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<",
                                           "db_file": ">>db_file<<",
                                           "spec": {"_priority": priority}}}
            fw2 = OptimizeFW(structure=structure, name=descriptor + 'rescale_optimize',
                             insert_db=False, job_type="normal", parents=[fw1],
                             **optimize_args["run_specs"], max_force_threshold=None)
            fw2 = powerups.add_cont_structure(fw2)
            fw2 = powerups.add_pass_structure(fw2, rescale_volume=True)

            if len(spawner_args["md_params"].keys()) > 0:
                run_args["md_params"].update(spawner_args["md_params"])
            fw3 = MDFW(structure=structure, name=descriptor + 'run1', previous_structure=True,
                       insert_db=False, **run_args["md_params"], parents=[fw2],
                       **run_args["run_specs"], **run_args["optional_fw_params"])

            fw3 = powerups.add_converge_task(fw3, **_spawner_args)

            fw_list.extend([fw1, fw2, fw3])
        else:
            fw1 = MDFW(structure=structure, name=descriptor + 'run0', previous_structure=False,
                       insert_db=False, parents=parents, **run_args["md_params"],
                       **run_args["run_specs"], **run_args["optional_fw_params"])
            fw1 = powerups.add_converge_task(fw1, **_spawner_args)
            fw_list.extend([fw1])

    # Production length MD runs
    # TODO Build continuation of MD on FIZZLED(from walltime) firework
    prod_steps = 0
    i = 0
    tag_id = uuid.uuid4()
    while prod_steps <= target_steps - max_steps:
        run_args = {"md_params": {"start_temp": run_args["md_params"]["end_temp"],
                                  "end_temp": run_args["md_params"]["end_temp"],
                                  "nsteps": max_steps},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<",
                                  "db_file": ">>db_file<<", "wall_time": 86400},
                    "optional_fw_params": {"override_default_vasp_params": {},
                                           "copy_vasp_outputs": False, "spec": {}},
                    "label": "prod_run"
                    }

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'PREC': 'Normal'}})
        run_args = recursive_update(run_args, prod_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority
        parents = fw_list[-1] if len(fw_list) > 0 else parents
        previous_structure = False if preconverged and i == 0 else True

        name = '-'.join([run_args["label"], str(i), str(tag_id)])
        fw = MDFW(structure=structure, name=name, previous_structure=previous_structure,
                  insert_db=True, parents=parents, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"])
        fw_list.append(fw)
        prod_steps += max_steps
        i += 1

    return fw_list


def get_converge_new(structure, temperature, converge_scheme='EOS', priority=None,
                     max_steps=5000, target_steps=10000, preconverged=False,
                     parents=None, **kwargs):
    """

    :param structure: Starting structure for the run
    :param temperature: Temperature for which to obtain a liquid
    :param images: Perturbations to the volume for the fit (fraction of input structure volume)
    :param preconverged: Is the structure already converged (i.e. Pressure 0bar) or volume rescaling not desired?
    :param max_steps: Maximum number of steps per chunk of production run MD simulation
    :param target_steps: Target number of steps for production MD run
    :param convergence_criteria: Type of convergence can specify ("density", <tolerance value>) or ("ionic", <tolerance value>)
    :param priority: Priority of all fireworks in the workflow
    :param kwargs: arguments such as spawner_args, converge_args, convergence_criteria, etc.
    :return:
    """
    # Generate a unique identifier for the fireworks belonging to this workflow
    tag_id = uuid.uuid4()

    fw_list = []

    # Setup initial Run and convergence of structure
    run_args = {"md_params": {"start_temp": temperature, "end_temp": temperature, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<",
                              "db_file": ">>db_file<<", "wall_time": 86400},
                "optional_fw_params": {
                    "override_default_vasp_params":
                        {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}},
                    "spec": {}}}
    run_args["optional_fw_params"]["spec"]["_queueadapter"] = {"walltime": run_args["run_specs"]["wall_time"]}
    run_args = recursive_update(run_args, kwargs.get('converge_args', {}))
    run_args["optional_fw_params"]["spec"]["_priority"] = priority

    # Setup Dictionary specifying parameters of the spawner for convergence tasks
    _spawner_args = {
        "converge_params": {"max_rescales": 15, "spawn_count": 1},
        "rescale_params": {"beta": 5e-7},
        "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
        "optional_fw_params": run_args["optional_fw_params"],
        "unique_identifier": tag_id}
    _spawner_args['converge_params']['converge_type'] = kwargs.get('convergence_criteria',
                                                                   [("density", 5), ('ionic', 0.001)])
    _spawner_args["md_params"].update({"start_temp": run_args["md_params"]["end_temp"]})
    _spawner_args["optional_fw_params"]["spec"]["_priority"] = priority
    _spawner_args = recursive_update(_spawner_args, kwargs.get('spawner_args', {}))

    # Converge the pressure (volume) of the system
    if not preconverged:
        if converge_scheme == 'EOS':
            # Create structures for varying volumes
            images = kwargs.get('images', [0.8, 1, 1.2])
            structures = [structure.copy() for i in images]
            for i, factor in enumerate(images):
                structures[i].scale_lattice(structure.volume * factor)

            # Create firework for each structure
            EOS_run_args = deepcopy(run_args)
            EOS_run_args['md_params']['nsteps'] = 1500 # Use less steps... Pressure usually converges rapidly
            EOS_run_args = recursive_update(EOS_run_args, kwargs.get('converge_args', {}))
            volume_fws = []
            for n, (i, vol_structure) in enumerate(zip(images, structures)):
                save_structure = True if n == len(images) - 1 else False
                _fw = MDFW(structure=vol_structure, name="volume_" + str(i) + "-" + str(tag_id),
                           previous_structure=False, insert_db=False,
                           **EOS_run_args["md_params"], **EOS_run_args["run_specs"],
                           **EOS_run_args["optional_fw_params"])

                _fw = powerups.add_pass_pv(_fw)
                volume_fws.append(_fw)
            fw_list.extend(volume_fws)

            # Create firework to converge pressure/volume
            _spawner_args['rescale_params']['beta'] = 1e-6
            spawner_fw = MDFW(structure=structure, name="run1" + "-" + str(tag_id),
                              previous_structure=True, insert_db=False,
                              parents=volume_fws, **run_args["md_params"],
                              **run_args["run_specs"],
                              **run_args["optional_fw_params"])

            spawner_fw = powerups.add_pv_volume_rescale(spawner_fw)
            spawner_fw = powerups.add_pass_pv(spawner_fw)
            spawner_fw = powerups.add_converge_task(spawner_fw, **_spawner_args)
            fw_list.append(spawner_fw)
        elif converge_scheme == 'vasp_opt':
            fw1 = MDFW(structure=structure, name="run0" + "-" + str(tag_id), previous_structure=False, insert_db=False,
                       **run_args["md_params"], **run_args["run_specs"],
                       **run_args["optional_fw_params"])

            # Create Vasp Optimization firework
            optimize_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                           "spec": {"_priority": priority}}}
            fw2 = OptimizeFW(structure=structure, name="rescale_optimize" + "-" + str(tag_id), insert_db=False,
                             job_type="normal",
                             parents=[fw1], **optimize_args["run_specs"], max_force_threshold=None)
            fw2 = powerups.add_cont_structure(fw2)
            fw2 = powerups.replace_pass_structure(fw2, rescale_volume=True)

            # Create Convergence Firework
            if len(_spawner_args["md_params"].keys()) > 0:
                run_args["md_params"].update(_spawner_args["md_params"])
            fw3 = MDFW(structure=structure, name="run1" + "-" + str(tag_id), previous_structure=True, insert_db=False,
                       **run_args["md_params"],
                       parents=[fw2], **run_args["run_specs"], **run_args["optional_fw_params"])

            fw3 = powerups.add_converge_task(fw3, **_spawner_args)

            fw_list.extend([fw1, fw2, fw3])
        else:
            fw1 = MDFW(structure=structure, name="run0" + "-" + str(tag_id), previous_structure=False, insert_db=False,
                       **run_args["md_params"], **run_args["run_specs"],
                       **run_args["optional_fw_params"])
            fw1 = powerups.add_converge_task(fw1, **_spawner_args)
            fw_list.append(fw1)

    # Production length MD runs
    prod_steps = 0
    i = 0
    while prod_steps <= target_steps - max_steps:
        run_args = {"md_params": {"start_temp": run_args["md_params"]["end_temp"],
                                  "end_temp": run_args["md_params"]["end_temp"],
                                  "nsteps": max_steps},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<",
                                  "db_file": ">>db_file<<", "wall_time": 86400},
                    "optional_fw_params": {"override_default_vasp_params": {},
                                           "copy_vasp_outputs": False, "spec": {}},
                    "label": "prod_run"
                    }

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'PREC': 'Normal'}})
        run_args["optional_fw_params"]["spec"]["_priority"] = priority
        parents = fw_list[-1] if len(fw_list) > 0 else parents
        previous_structure = False if preconverged and i == 0 else True

        name = '-'.join([run_args["label"], str(i), str(tag_id)])
        fw = MDFW(structure=structure, name=name, previous_structure=previous_structure,
                  insert_db=True, parents=parents, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"])
        fw_list.append(fw)
        prod_steps += max_steps
        i += 1

    return fw_list
