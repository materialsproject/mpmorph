from fireworks import Firework
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.fireworks.core import StaticFW
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from mpmorph.firetasks.mdtasks import ConvergeTask
from mpmorph.firetasks.glue_tasks import PreviousStructureTask, SaveStructureTask
from mpmorph.sets import FrozenPhononSet
from mpmorph.firetasks.dbtasks import VaspMDToDb
from pymatgen.io.vasp.sets import MITMDSet, MPStaticSet


class MDFW(Firework):
    def __init__(self, structure, start_temp, end_temp, nsteps, name="molecular dynamics",
                 vasp_input_set=None, vasp_cmd="vasp", override_default_vasp_params=None,
                 wall_time=19200, db_file=None, parents=None, copy_vasp_outputs=True,
                 previous_structure=False, insert_db=False, **kwargs):
        """
        This Firework is modified from atomate.vasp.fireworks.core.MDFW to fit the needs of mpmorph
        Standard firework for a single MD run.
        Args:
            structure (Structure): Input structure.
            start_temp (float): Start temperature of MD run.
            end_temp (float): End temperature of MD run.
            nsteps (int): Number of MD steps
            name (string): Name for the Firework.
            vasp_input_set (string): string name for the VASP input set (e.g.,
                "MITMDVaspInputSet").
            vasp_cmd (string): Command to run vasp.
            override_default_vasp_params (dict): If this is not None,
                these params are passed to the default vasp_input_set, i.e.,
                MITMDSet. This allows one to easily override some
                settings, e.g., user_incar_settings, etc. Particular to MD,
                one can control time_step and all other settings of the input set.
            wall_time (int): Total wall time in seconds before writing STOPCAR.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MITMDSet(structure, start_temp=start_temp,
                                                    end_temp=end_temp, nsteps=nsteps,
                                                    **override_default_vasp_params)

        t = []
        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                     contcar_to_poscar=True))

        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        if previous_structure:
            t.append(PreviousStructureTask())
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time))
        t.append(PassCalcLocs(name=name))
        t.append(SaveStructureTask())
        if insert_db:
            t.append(VaspMDToDb(db_file=db_file,
                              additional_fields={"task_label": name}, defuse_unsuccessful=False))
        super(MDFW, self).__init__(t, parents=parents,
                                   name="{}-{}".format(structure.composition.reduced_formula, name),
                                   **kwargs)

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

class StaticFW(Firework):
    def __init__(self, structure, name="static", vasp_input_set=None, vasp_cmd="vasp",
                 db_file=None, parents=None, **kwargs):
        """
        This Firework is modified from atomate.vasp.fireworks.core.StaticFW to fit the needs of mpmorph
        Standard static calculation Firework - either from a previous location or from a structure.
        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, grabs a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        t = []
        vasp_input_set = vasp_input_set or MPStaticSet(structure)
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)