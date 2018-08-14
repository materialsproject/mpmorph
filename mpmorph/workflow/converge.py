from fireworks import Workflow
from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW, OptimizeFW
from mpmorph.util import recursive_update
import uuid


def get_converge(structure, priority=None, preconverged=False, max_steps=5000, target_steps=40000,
                 spawner_args={}, converge_args={}, prod_args={}, converge_type=("density", 5), vasp_opt=True, **kwargs):
    """

    :param structure:
    :param temperatures:
    :param priority:
    :param preconverged: Is the structure already converged (i.e. Pressure 0bar) or volume rescaling not desired?
    :param prod_quants:
    :param spawner_args:
    :param converge_args:
    :param prod_args:
    :param converge_type:
    :param vasp_opt:
    :param kwargs:
    :return:
    """
    fw_list = []
    # Initial Run and convergence of structure

    run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "wall_time": 86400},
                "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
    run_args["optional_fw_params"]["spec"]["_queueadapter"] = {"walltime": run_args["run_specs"]["wall_time"]}
    run_args = recursive_update(run_args, converge_args)
    run_args["optional_fw_params"]["spec"]["_priority"] = priority

    _spawner_args = {"converge_params": {"converge_type": [converge_type], "max_rescales": 15, "spawn_count": 1},
                     "rescale_params": {"beta": 0.0000005},
                     "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
                     "optional_fw_params": run_args["optional_fw_params"]}
    _spawner_args["md_params"].update({"start_temp": run_args["md_params"]["end_temp"]})
    _spawner_args = recursive_update(_spawner_args, spawner_args)

    if not preconverged:
        if vasp_opt:
            fw1 = MDFW(structure=structure, name="run0", previous_structure=False, insert_db=False,
                       **run_args["md_params"], **run_args["run_specs"],
                       **run_args["optional_fw_params"])

            # OptimizeFW does not take wall_time
            optimize_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                           "spec": {"_priority": priority}}}
            fw2 = OptimizeFW(structure=structure, name="rescale_optimize", insert_db=False, job_type="normal",
                             parents=[fw1], **optimize_args["run_specs"], max_force_threshold=None)
            fw2 = powerups.add_cont_structure(fw2)
            fw2 = powerups.replace_pass_structure(fw2, rescale_volume=True)

            if len(spawner_args["md_params"].keys()) > 0:
                run_args["md_params"].update(spawner_args["md_params"])
            fw3 = MDFW(structure=structure, name="run1", previous_structure=True, insert_db=False, **run_args["md_params"],
                       parents=[fw2], **run_args["run_specs"], **run_args["optional_fw_params"])

            fw3 = powerups.add_converge_task(fw3, **_spawner_args)

            fw_list.extend([fw1, fw2, fw3])
        else:
            fw1 = MDFW(structure=structure, name="run0", previous_structure=False, insert_db=False,
                       **run_args["md_params"], **run_args["run_specs"],
                       **run_args["optional_fw_params"])
            fw1 = powerups.add_converge_task(fw1, **_spawner_args)
            fw_list.extend([fw1])

    # Production length MD runs
    # TODO Build continuation of MD on FIZZLED(from walltime) firework
    prod_steps = 0
    i = 0
    tag_id = uuid.uuid4()
    while prod_steps <= target_steps - max_steps:
        run_args = {"md_params": {"start_temp": run_args["md_params"]["end_temp"],
                                  "end_temp": run_args["md_params"]["end_temp"], "nsteps": max_steps},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                  "wall_time": 86400},
                    "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}},
                    "label": "prod_run_"}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args = recursive_update(run_args, prod_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority
        parents = fw_list[-1] if len(fw_list) > 0 else []
        previous_structure = False if preconverged and i == 0 else True
        fw = MDFW(structure=structure, name=run_args["label"] + str(i) + "-" + str(tag_id),
                  previous_structure=previous_structure, insert_db=True, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"], parents=parents)
        fw_list.append(fw)

        prod_steps += max_steps
        i += 1

    # if to_trajectory:
    #     fw_list[-1] = powerups.aggregate_trajectory(fw_list[-1], tag_id)

    pretty_name = structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name=pretty_name + "_diffusion")
    return wf
