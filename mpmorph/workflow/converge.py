from fireworks import Firework, Workflow
from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW
from mpmorph.util import recursive_update

def get_converge(structure, priority = None, preconverged=False, prod_quants={"nsteps":5000,"target": 40000}, spawner_args={}, converge_args={}, prod_args={}, converge_type=("density", 5), **kwargs):
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
    :param kwargs:
    :return:
    """
    fw_list = []
    #Initial Run and convergence of structure

    run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "wall_time": 86400},
                "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}}}
    if not preconverged:

        run_args["optional_fw_params"]["override_default_vasp_params"].update({'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args = recursive_update(run_args, converge_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority

        fw = MDFW(structure=structure, name = "run0", previous_structure=False, insert_db=False, **run_args["md_params"],**run_args["run_specs"], **run_args["optional_fw_params"])

        _spawner_args = {"converge_params":{"converge_type": [converge_type], "max_rescales": 10, "spawn_count": 0},
                         "rescale_args":{"beta": 0.0000008},
                         "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
                         "optional_fw_params":run_args["optional_fw_params"]}
        _spawner_args["md_params"].update({"start_temp":run_args["md_params"]["end_temp"]})
        _spawner_args = recursive_update(_spawner_args, spawner_args)

        fw = powerups.add_converge_task(fw, **_spawner_args)

        fw_list.append(fw)

    #Production length MD runs
    #TODO Build continuation of MD on FIZZLED(from walltime) firework
    prod_steps = 0
    i = 0
    while prod_steps <= prod_quants["target"] - prod_quants["nsteps"]:
        run_args = {"md_params": {"start_temp": run_args["md_params"]["end_temp"], "end_temp": run_args["md_params"]["end_temp"], "nsteps":5000},
                    "run_specs":{"vasp_input_set": None ,"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<", "wall_time": 86400},
                    "optional_fw_params":{"override_default_vasp_params":{}, "copy_vasp_outputs": False, "spec":{}},
                    "label": "prod_run_"}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args = recursive_update(run_args, prod_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority
        parents = fw_list[-1] if len(fw_list) > 0 else []
        fw = MDFW(structure=structure, name = run_args["label"] + str(i), previous_structure=True, insert_db=True, **run_args["md_params"], **run_args["run_specs"], **run_args["optional_fw_params"], parents=parents)
        fw_list.append(fw)

        prod_steps += prod_quants["nsteps"]
        i+=1

    pretty_name=structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name = pretty_name + "_diffusion")
    return wf
