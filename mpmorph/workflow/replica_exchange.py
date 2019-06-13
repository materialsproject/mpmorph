from fireworks import Firework, Workflow
from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW, OptimizeFW
from mpmorph.util import recursive_update
import uuid


def get_temper(structures, swap_frequency=500, priority=None):
    """

    :param structure: Starting structure for all MD simulations
    :param temperatures: List of temperatures for which replica exchange is performed
    :param kwargs:
    :return:
    """

    temperatures = []
    fw_list = []

    run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "wall_time": 86400},
                "optional_fw_params": {
                    "override_default_vasp_params": {{'user_incar_settings': {'ISIF': 1, 'LWAVE': False}}},
                    "copy_vasp_outputs": False,
                    "spec": {"_priority": priority, "_queueadapter": {"walltime": 86400}}}}


    # TODO Setup parallel tempering run

    #

    pretty_name = structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name=pretty_name + "_temper")
    return wf
