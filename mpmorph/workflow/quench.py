from fireworks import Firework, Workflow
from pymatgen import Structure, Composition
from mpmorph.fireworks import powerups
from atomate.vasp.fireworks.core import MDFW, StaticFW, OptimizeFW
from mpmorph.util import recursive_update
import numpy as np


def get_quench(structures, temperatures={}, priority = None, type="simulated_anneal", cool_args = {}, hold_args = {}, unconverged=True, **kwargs):
    fw_list = []
    if temperatures == {}:
        temperatures = {"start_temp": 3000, "end_temp": 500, "temp_step":500}
    if cool_args == {}:
        cool_args = {"md_params": {"nsteps": 200}}
    if hold_args == {}:
        hold_args = {"md_params": {"nsteps": 500}}

    for (i, structure) in enumerate(structures):
        _fw_list=[]
        if type == "simulated_anneal":
            for temp in np.arange(temperatures["start_temp"], temperatures["end_temp"]-temperatures["temp_step"], -temperatures["temp_step"]):
                # get fw for cool step
                _fw = get_MDFW(structure, temp, temp-temperatures["temp_step"], args=cool_args,
                                         parents=_fw_list[-1] if len(_fw_list) > 0 else [])
                _fw = powerups.replace_vaspmdtodb(_fw)
                _fw = powerups.add_pass_structure(_fw)
                if len(_fw_list) > 0:
                    _fw = powerups.add_cont_structure(_fw)
                _fw_list.append(_fw)
                # get fw for hold step
                _fw =get_MDFW(structure, temp-temperatures["temp_step"], temp-temperatures["temp_step"],
                                         args=hold_args, parents=_fw_list[-1])
                _fw = powerups.replace_vaspmdtodb(_fw)
                _fw = powerups.add_pass_structure(_fw)
                _fw = powerups.add_cont_structure(_fw)
                _fw_list.append(_fw)

            fw_list.extend(_fw_list)

        # Relax OptimizeFW and StaticFW
        run_args = {"run_specs":{"vasp_input_set": None ,"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<", "wall_time": 86400, "spec":{}},
                    "optional_fw_params":{"override_default_vasp_params":{}, "copy_vasp_outputs": False}}
        _name = str(structure.composition.reduced_formula)+"_"+str(i)

        fw1 = OptimizeFW(structure=structure, name=_name + "optimize", parents= _fw_list[-1] if len(_fw_list) > 0 else [], **run_args["run_specs"])
        if len(_fw_list) > 0:
            fw1 = powerups.add_cont_structure(fw1)
        fw1 = powerups.add_pass_structure(fw1)

        fw2 = StaticFW(structure=structure, name=_name + "static", parents=[fw1], **run_args["run_specs"])
        fw2 = powerups.add_cont_structure(fw2)
        fw2 = powerups.add_pass_structure(fw2)

        _fw_list.extend([fw1, fw2])
    wf = Workflow(fw_list, "")

def get_MDFW(structure, start_temp, end_temp, priority=None, args={}, **kwargs):
    run_args = {"md_params": {"nsteps": 500},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "wall_time": 40000},
                "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}}}

    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
    run_args = recursive_update(run_args, args)
    run_args["md_params"]["start_temp"] = start_temp
    run_args["md_params"]["end_temp"] = end_temp
    run_args["optional_fw_params"]["spec"]["_priority"] = priority
    _mdfw = MDFW(structure=structure, name = "run0", **run_args["md_params"],**run_args["run_specs"], **run_args["optional_fw_params"])
    return _mdfw