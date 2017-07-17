from fireworks import Firework
from atomate.vasp.fireworks.core import MDFW

class ConvergeFW(Firework):
    """
    Atomate MDFW with additional SpawnMDFW task added to end
    """
    def __init__(self):
        mdfw = MDFW()
        spawner_task = SpawnMDFW()
        mdfw.tasks.append(spawner_task)
        return mdfw

class DiffusionFW(Firework):
    pass

class QuenchFW(Firework):
    pass
