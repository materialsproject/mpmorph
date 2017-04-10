from fireworks import explicit_serialize, FireTaskBase, FWAction, Firework, LaunchPad
from mpmorph.runners.amorphous_maker import AmorphousMaker
from mpmorph.runners.rescale_volume import RescaleVolume
from mpmorph.analysis.md_data import parse_pressure
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet, ModifyIncar
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MITMDSet
from atomate.common.firetasks.glue_tasks import PassCalcLocs
import shutil
import numpy as np
from atomate.vasp.firetasks.parse_outputs import VaspToDbTask
import os

__author__ = 'Muratahan Aykol <maykol@lbl.gov>'


# TODO: 2. Add option to lead to a production run of specified length after density is found
# TODO: 3. Switch to MPRelax Parameters in MD
# TODO: 4. Database insertion?
# TODO: 5. Parser tasks

@explicit_serialize
class AmorphousMakerTask(FireTaskBase):
    """
    Create a constrained-random packed structure from composition and box dimensions.
    Required params:
        composition: (dict) a dict of target composition with integer atom numbers
                        e.g. {"V":22, "Li":10, "O":75, "B":10}
        box_scale: (float) all lattice vectors are multiplied with this scalar value.
                        e.g. edge length of a cubic simulation box.
    Optional params:
        tol (float): tolerance factor for how close the atoms can get (angstroms).
                        e.g. tol = 2.0 angstroms
        packmol_path (str): path to the packmol executable. Defaults to "packmol"
        clean (bool): whether the intermedite files generated are deleted. Defaults to True.
    """

    required_params = ["composition", "box_scale"]
    optional_params = ["packmol_path", "clean", "tol"]

    def run_task(self, fw_spec):
        glass = AmorphousMaker(self.get("composition"), self.get("box_scale"), self.get("tol", 2.0),
                               packmol_path=self.get("packmol_path", "packmol"),
                               clean=self.get("clean", True))
        structure = glass.random_packed_structure.as_dict()
        return FWAction(stored_data=structure)


@explicit_serialize
class GetPressureTask(FireTaskBase):
    required_params = ["outcar_path"]
    optional_params = ["averaging_fraction"]

    def run_task(self, fw_spec):
        p = parse_pressure(self["outcar_path"], self.get("averaging_fraction", 0.5))
        if fw_spec['avg_pres']:
            fw_spec['avg_pres'].append(p[0] * 1000)
        else:
            fw_spec['avg_pres'] = [p[0] * 1000]
        return FWAction()


@explicit_serialize
class SpawnMDFWTask(FireTaskBase):
    """
    Decides if a new MD calculation should be spawned or if density is found. If so, spawns a new calculation.
    """
    required_params = ["pressure_threshold", "max_rescales", "vasp_cmd", "wall_time",
                       "db_file", "spawn_count", "copy_calcs", "calc_home"]
    optional_params = ["averaging_fraction", "cool"]

    def run_task(self, fw_spec):
        vasp_cmd = self["vasp_cmd"]
        wall_time = self["wall_time"]
        db_file = self["db_file"]
        max_rescales = self["max_rescales"]
        pressure_threshold = self["pressure_threshold"]
        spawn_count = self["spawn_count"]
        calc_home = self["calc_home"]
        copy_calcs = self["copy_calcs"]
        snaps = self["cool"] or False

        if spawn_count > max_rescales:
            # TODO: Log max rescale reached info.
            _temp_var = 1
            #return FWAction(defuse_workflow=True)

        name = ("spawnrun" + str(spawn_count))

        current_dir = os.getcwd()
        averaging_fraction = self.get("averaging_fraction", 0.5)
        p = parse_pressure("./", averaging_fraction)[0]

        if np.fabs(p) > pressure_threshold:
            t = []
            # Copy the VASP outputs from previos run. Very first run get its from the initial MDWF which
            # uses PassCalcLocs. For the rest we just specify the previous dir.
            if spawn_count == 0:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=False))
            else:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=True))

            t.append(RescaleVolumeTask(initial_pressure=p * 1000.0, initial_temperature=1))
            t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                      handler_group="md", wall_time=wall_time, gzip_output=False))
            # Will implement the database insertion
            # t.append(VaspToDbTask(db_file=db_file,
            #                       additional_fields={"task_label": "density_adjustment"}))
            if copy_calcs:
                t.append(CopyCalsHome(calc_home=calc_home, run_name=name))
            t.append(SpawnMDFWTask(pressure_threshold=pressure_threshold,
                                   max_rescales=max_rescales,
                                   wall_time=wall_time,
                                   vasp_cmd=vasp_cmd,
                                   db_file=db_file,
                                   spawn_count=spawn_count + 1,
                                   copy_calcs=copy_calcs,
                                   calc_home=calc_home,
                                   averaging_fraction=averaging_fraction,
                                   cool=snaps))
            new_fw = Firework(t, name=name)
            return FWAction(stored_data={'pressure': p}, additions=[new_fw])

        else:
            t = []
            name = "long_melt"
            if spawn_count == 0:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=False))
            else:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=True))

            t.append(ModifyIncar({'NSW': 10000}))
            t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time, gzip_output=False))

            if copy_calcs:
                t.append(CopyCalsHome(calc_home=calc_home, run_name=name))

            t.append(PassCalcLocs(name=name))
            if snaps:
                t.append(StructureSamplerTask(copy_calcs=copy_calcs, calc_home=calc_home))
                fw = Firework(t, name=name)
                return FWAction(stored_data={'pressure': p, 'density_calculated': True}, additions={fw})
            fw = Firework(t, name=name)
            return FWAction(additions={fw})


@explicit_serialize
class RescaleVolumeTask(FireTaskBase):
    """
    Volume rescaling
    """
    required_params = ["initial_temperature", "initial_pressure"]
    optional_params = ["target_pressure", "target_temperature", "target_pressure", "alpha", "beta"]

    def run_task(self, fw_spec):
        # Initialize volume correction object with last structure from last_run
        initial_temperature = self["initial_temperature"]
        initial_pressure = self["initial_pressure"]
        target_temperature = self.get("target_temperature", initial_temperature)
        target_pressure = self.get("target_pressure", 0.0)
        alpha = self.get("alpha", 10e-6)
        beta = self.get("beta", 10e-7)
        corr_vol = RescaleVolume.of_poscar(poscar_path="./POSCAR", initial_temperature=initial_temperature,
                                           initial_pressure=initial_pressure,
                                           target_pressure=target_pressure,
                                           target_temperature=target_temperature, alpha=alpha, beta=beta)
        # Rescale volume based on temperature difference first. Const T will return no volume change:
        corr_vol.by_thermo(scale='temperature')
        # TO DB ("Rescaled volume due to delta T: ", corr_vol.structure.volume)
        # Rescale volume based on pressure difference:
        corr_vol.by_thermo(scale='pressure')
        # TO DB ("Rescaled volume due to delta P: ", corr_vol.structure.volume)
        corr_vol.poscar.write_file("./POSCAR")
        # Pass the rescaled volume to Poscar
        return FWAction(stored_data=corr_vol.structure.as_dict())


@explicit_serialize
class CopyCalsHome(FireTaskBase):
    required_params = ["calc_home", "run_name"]
    optional_params = ["files"]

    def run_task(self, fw_spec):
        default_list = ["INCAR", "POSCAR", "CONTCAR", "OUTCAR", "POTCAR", "vasprun.xml", "XDATCAR", "OSZICAR", "DOSCAR"]
        files = self.get("files", default_list)
        calc_home = self["calc_home"]
        run_name = self["run_name"]
        target_dir = os.path.join(calc_home, run_name)
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)
        for f in files:
            try:
                shutil.copy2(f, target_dir)
            except:
                pass
        return FWAction()


@explicit_serialize
class SimulatedAnnealTask(FireTaskBase):
    """
    Spawn fireworks until temp threshold
    """
    required_params = ["temperature", "temperature_threshold", "temperature_decrement", "vasp_cmd", "wall_time",
                       "db_file", "copy_calcs", "calc_home", "spawn_count", "nsteps_cool", "nsteps_hold", "name"]
    optional_params = []

    def run_task(self, fw_spec):
        temperature = self["temperature"]
        temperature_threshold = self["temperature_threshold"]
        temperature_decrement = self["temperature_decrement"]
        vasp_cmd = self["vasp_cmd"]
        wall_time = self["wall_time"]
        db_file = self["db_file"]
        copy_calcs = self["copy_calcs"]
        calc_home = self["calc_home"]
        spawn_count = self["spawn_count"]
        name = self["name"]

        # Check Temperature
        if temperature >= temperature_threshold:
            fw1 = self.CoolMD()
            t=SimulatedAnnealTask(temperature=temperature - temperature_decrement,
                                         temperature_threshold=temperature_threshold,
                                         temperature_decrement=temperature_decrement,
                                         vasp_cmd=vasp_cmd,
                                         wall_time=wall_time,
                                         db_file=db_file,
                                         copy_calcs=copy_calcs,
                                         calc_home=calc_home,
                                         spawn_count=spawn_count + 1, )
            fw2 = self.HoldMD(sat=t)
            return FWAction(additions=[fw1, fw2])
        else:
            fws = []
            fw1 = OptimizeFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[], **kwargs)
            fws.append(fw1)
            fw2 = StaticFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[fw1])
            fws.append(fw2)
            if copy_calcs:
                t=[CopyCalsHome(os.path.join(calc_home, name), run_name="final_relax")]
                fw3=Firework(t, name=name+"_copy_calcs")
            return FWAction()

    def CoolMD(self, nsteps=None):
        temperature = self["temperature"]
        temperature_decrement = self["temperature_decrement"]
        vasp_cmd = self["vasp_cmd"]
        wall_time = self["wall_time"]
        copy_calcs = self["copy_calcs"]
        calc_home = self["calc_home"]
        nsteps_cool = nsteps or self["nsteps_cool"]
        name = self["name"]

        t = []
        t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))
        t.append(ModifyIncar({'TEBEG': temperature, 'TEEND': temperature - temperature_decrement, 'NSW': nsteps_cool}))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time))
        t.append(PassCalcLocs(name=name+"_cool_"+str(temperature-temperature_decrement)))
        if copy_calcs:
            t.append(CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name=name+"_cool_"+str(temperature-temperature_decrement)))

        return Firework(t, name=name+"_cool_"+str(temperature-temperature_decrement))

    def HoldMD(self, nsteps=None, sat=None):
        temperature = self["temperature"]
        temperature_decrement = self["temperature_decrement"]
        vasp_cmd = self["vasp_cmd"]
        wall_time = self["wall_time"]
        copy_calcs = self["copy_calcs"]
        calc_home = self["calc_home"]
        nsteps_heat = nsteps or self["nsteps_heat"]
        name = self["name"]

        t = []
        t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))
        t.append(ModifyIncar({'TEBEG': temperature - temperature_decrement, 'TEEND': temperature - temperature_decrement,
             'NSW': nsteps_heat}))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time))
        t.append(PassCalcLocs(name=name+"_hold_"+str(temperature-temperature_decrement)))
        if copy_calcs:
            t.append(CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name=name+"_hold_"+str(temperature-temperature_decrement)))
        t.append(sat)
        return Firework(t, name=name+"_hold_"+str(temperature-temperature_decrement))

@explicit_serialize
class StructureSamplerTask(FireTaskBase):
    """

    """
    required_params = ["copy_calcs", "calc_home"]
    optional_params = []

    def run_task(self, fw_spec):
        copy_calcs = self["copy_calcs"]
        calc_home = self["calc_home"]
        current_dir = os.getcwd()
        xdatcar_file = os.path.join(current_dir, 'XDATCAR')
        wfs = get_wf_structure_sampler(xdatcar_file=xdatcar_file, sim_anneal=True, copy_calcs=copy_calcs, calc_home=calc_home, n=1)
        lp = LaunchPad()
        for _wf in wfs:
            lp.add_wf(_wf)
        return FWAction()


@explicit_serialize
class VaspMdToDbTask(FireTaskBase):
    pass


@explicit_serialize
class VaspMdToDiffusion(FireTaskBase):
    pass


@explicit_serialize
class VaspMdToStructuralAnalysis(FireTaskBase):
    pass

from mpmorph.workflow.workflows import get_wf_structure_sampler