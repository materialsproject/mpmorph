from fireworks import Firework, Workflow
from mpmorph.fireworks import powerups
from atomate.vasp.fireworks.core import MDFW


def get_converge(structure, priority = None, preconverged=False, prod_quants={"nsteps":5000,"target": 40000}, spawner_args={}, converge_args={}, prod_args={}, converge_type="density", **kwargs):
    """

    :param structure:
    :param temperatures:
    :param priority:
    :param preconverged:
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

    if not preconverged:
        run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps":2000},
                    "run_specs":{"vasp_input_set": None ,"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<", "wall_time": 86400},
                    "optional_fw_params":{"override_default_vasp_params":{}, "copy_vasp_outputs": False, "spec":{}}}

        run_args["optional_fw_params"]["override_default_vasp_params"].update({'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args.update(converge_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority

        fw = MDFW(structure=structure, name = "run0", **run_args["md_params"],**run_args["run_specs"], **run_args["optional_fw_params"])

        _spawner_args = {"converge_params":{"converge_type": [("density", 5)], "max_rescales": 10, "spawn_count": 0},
                         "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
                         "optional_fw_params":run_args["optional_fw_params"]}
        _spawner_args["md_params"].update({"start_temp":run_args["md_params"]["end_temp"]})
        _spawner_args.update(_spawner_args)

        fw = powerups.replace_vaspmdtodb(fw)
        fw = powerups.add_converge_task(fw, **_spawner_args)
        fw = powerups.add_pass_structure(fw)

        fw_list.append(fw)

    #Production length MD runs
    #TODO Build continuation of MD on FIZZLED(from walltime) firework
    prod_steps = 0
    i = 0
    while prod_steps <= prod_quants["target"] - prod_quants["nsteps"]:
        run_args = run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps":2000},
                    "run_specs":{"vasp_input_set": None ,"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<", "wall_time": 86400},
                    "optional_fw_params":{"override_default_vasp_params":{}, "copy_vasp_outputs": False, "spec":{}}}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args.update(prod_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority
        parents = fw_list[-1] if len(fw_list) > 0 else []
        fw = MDFW(structure=structure, name = "prod_run_" + str(i), **run_args["md_params"], **run_args["run_specs"], **run_args["optional_fw_params"], parents=parents)
        fw = powerups.replace_vaspmdtodb(fw)
        fw = powerups.add_cont_structure(fw, position=1) #Add after MDFW WriteInputSet to override structure
        fw = powerups.add_pass_structure(fw) #Add at end to append structure to fw_spec
        fw_list.append(fw)

        prod_steps += prod_quants["nsteps"]
        i+=1

    pretty_name=structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name = pretty_name + "_diffusion")
    return wf