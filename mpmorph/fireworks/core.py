from fireworks import Firework
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.fireworks.core import StaticFW
from mpmorph.firetasks.mdtasks import ConvergeTask
from mpmorph.sets import FrozenPhononSet

from atomate.vasp.firetasks.parse_outputs import VaspToDb

class ConvergeFW(Firework):
    pass

class DiffusionFW(Firework):
    pass

class QuenchFW(Firework):
    pass

class VibrationalFW(Firework):
    def __init__(self, structure, name="vibrational_fw", vasp_cmd=">>vasp_cmd<<",
                 db_file=">>db_file<<", incar_updates={}, parents=None, **kwargs):
        # Write input set for vibrational
        t=[]
        vasp_input_set = FrozenPhononSet(structure, incar_updates)
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        # Job tasks
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)