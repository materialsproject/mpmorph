from fireworks import Firework, Workflow
from pymatgen import Structure, Composition
from mpmorph.fireworks.core import ConvergeFW, DiffusionFW
import mpmorph.fireworks.powerups as mp_powerup
from atomate.vasp.fireworks.core import MDFW


def get_quench(structure, temperatures, type="simulated_anneal", converge_type="density", unconverged=True, **kwargs):
    fw_list = []

    if type == "simulated_anneal":
        #Do something
        fw = MDFW()
        mp_powerup.add_cont_structure(fw)
        fw_list.append(fw)
    else:
