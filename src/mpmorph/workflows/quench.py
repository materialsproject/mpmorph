import numpy as np
from fireworks import Workflow

from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW, OptimizeFW, StaticFW
from mpmorph.utils import recursive_update


def get_quench_wf(
    structures, priority=None, quench_type="slow_quench", descriptor="", **kwargs
):
    """

    Args:
        structure: Starting structure for the run
        priority: Priority of all fireworks in the workflows
        quench_type: use "slow_quench" for a gradual decrease in temperature or
            "mp_quench" for a instantaneous DFT relaxation
        target_steps: Target number of steps for production MD run
        descriptor: Extra description to add to the name of the firework
        **kwargs: Arguments such as cool_args, hold_args, quench_args, etc. Cool_args and hold args are only applicable
            when using "slow_quench"

    Returns: Workflow object

    """

    fw_list = []
    temperatures = kwargs.get(
        "temperatures", {"start_temp": 3000, "end_temp": 500, "temp_step": 500}
    )
    cool_args = kwargs.get("cool_args", {"md_params": {"nsteps": 200}})
    hold_args = kwargs.get("hold_args", {"md_params": {"nsteps": 500}})
    quench_args = kwargs.get("quench_args", {})

    for i, structure in enumerate(structures):
        _fw_list = []
        if quench_type == "slow_quench":
            for temp in np.arange(
                temperatures["start_temp"],
                temperatures["end_temp"],
                -temperatures["temp_step"],
            ):
                # get fw for cool step
                use_prev_structure = False
                if len(_fw_list) > 0:
                    use_prev_structure = True
                _fw = get_MDFW(
                    structure,
                    temp,
                    temp - temperatures["temp_step"],
                    name="snap_"
                    + str(i)
                    + "_cool_"
                    + str(temp - temperatures["temp_step"]),
                    args=cool_args,
                    parents=[_fw_list[-1]] if len(_fw_list) > 0 else [],
                    priority=priority,
                    previous_structure=use_prev_structure,
                    insert_db=True,
                    **kwargs,
                )
                _fw_list.append(_fw)
                # get fw for hold step
                _fw = get_MDFW(
                    structure,
                    temp - temperatures["temp_step"],
                    temp - temperatures["temp_step"],
                    name="snap_"
                    + str(i)
                    + "_hold_"
                    + str(temp - temperatures["temp_step"]),
                    args=hold_args,
                    parents=[_fw_list[-1]],
                    priority=priority,
                    previous_structure=True,
                    insert_db=True,
                    **kwargs,
                )
                _fw_list.append(_fw)

        if quench_type in ["slow_quench", "mp_quench"]:
            # Relax OptimizeFW and StaticFW
            run_args = {
                "run_specs": {
                    "vasp_input_set": None,
                    "vasp_cmd": ">>vasp_cmd<<",
                    "db_file": ">>db_file<<",
                    "spec": {"_priority": priority},
                },
                "optional_fw_params": {"override_default_vasp_params": {}},
            }
            run_args = recursive_update(run_args, quench_args)
            _name = "snap_" + str(i)

            use_prev_structure = True if len(_fw_list) > 0 else False
            fw1 = OptimizeFW(
                structure=structure,
                name=f"{_name}{descriptor}_optimize",
                parents=[_fw_list[-1]] if len(_fw_list) > 0 else [],
                previous_structure=use_prev_structure,
                **run_args["run_specs"],
                **run_args["optional_fw_params"],
                max_force_threshold=None,
            )

            fw2 = StaticFW(
                structure=structure,
                name=f"{_name}{descriptor}_static",
                parents=[fw1],
                previous_structure=True,
                **run_args["run_specs"],
                **run_args["optional_fw_params"],
            )

            _fw_list.extend([fw1, fw2])

        fw_list.extend(_fw_list)

    name = structure.composition.reduced_formula + descriptor + "_quench"
    wf = Workflow(fw_list, name=name)
    return wf


def get_MDFW(
    structure,
    start_temp,
    end_temp,
    name="molecular dynamics",
    priority=None,
    args={},
    **kwargs,
):
    """

    Helper function to get molecular dynamics firework for quench workflow

    Args:
        structure: Initial structure for molecular dynamics run
        start_temp: Starting Temperature
        end_temp: Ending Temperature
        name: name of firework
        priority: priority of job in database
        args: custom arguments dictionary for molecular dynamics run
        kwargs: kwargs for MDFW

    Returns: Molecular Dynamics Firework

    """
    # Get customized firework
    run_args = {
        "md_params": {"nsteps": 500, "start_temp": start_temp, "end_temp": end_temp},
        "run_specs": {
            "vasp_input_set": None,
            "vasp_cmd": ">>vasp_cmd<<",
            "db_file": ">>db_file<<",
            "wall_time": 40000,
        },
        "optional_fw_params": {
            "override_default_vasp_params": {},
            "spec": {"_priority": priority},
        },
    }

    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {"user_incar_settings": {"ISIF": 1, "LWAVE": False, "PREC": "Low"}}
    )
    run_args = recursive_update(run_args, args)
    _mdfw = MDFW(
        structure=structure,
        name=name,
        **run_args["md_params"],
        **run_args["run_specs"],
        **run_args["optional_fw_params"],
        **kwargs,
    )
    return _mdfw
