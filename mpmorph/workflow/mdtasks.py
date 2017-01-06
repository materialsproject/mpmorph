from fireworks import explicit_serialize, FireTaskBase, FWAction, Firework
from mpmorph.runners.amorphous_maker import AmorphousMaker
from mpmorph.runners.rescale_volume import RescaleVolume
from mpmorph.analysis.md_data import parse_pressure
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.common.firetasks.glue_tasks import PassCalcLocs
import shutil
import numpy as np
from atomate.vasp.firetasks.parse_outputs import VaspToDbTask
import os

__author__ = 'Muratahan Aykol <maykol@lbl.gov>'

#TODO: 2. Add option to lead to a production run of specified length after density is found
#TODO: 3. Switch to MPRelax Parameters in MD
#TODO: 4. Database insertion?
#TODO: 5. Parser tasks

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
            fw_spec['avg_pres'].append(p[0]*1000)
        else:
            fw_spec['avg_pres'] = [p[0]*1000]
        return FWAction()


@explicit_serialize
class SpawnMDFWTask(FireTaskBase):
    """
    Decides if a new MD calculation should be spawned or if density is found. If so, spawns a new calculation.
    """
    required_params = ["pressure_threshold", "max_rescales", "vasp_cmd", "wall_time",
                       "db_file", "spawn_count", "copy_calcs", "calc_home"]
    optional_params = ["averaging_fraction"]

    def run_task(self, fw_spec):
        vasp_cmd = self["vasp_cmd"]
        wall_time = self["wall_time"]
        db_file = self["db_file"]
        max_rescales = self["max_rescales"]
        pressure_threshold = self["pressure_threshold"]
        spawn_count = self["spawn_count"]
        calc_home = self["calc_home"]
        copy_calcs = self["copy_calcs"]

        if spawn_count > max_rescales:
            # TODO: Log max rescale reached info.
            return FWAction(defuse_workflow=True)

        name = ("spawnrun"+str(spawn_count))

        current_dir = os.getcwd()
        averaging_fraction = self.get("averaging_fraction", 0.5)
        p = parse_pressure("./", averaging_fraction)[0]

        if np.fabs(p) > pressure_threshold:
            t = []
            # Copy the VASP outputs from previos run. Very first run get its from the initial MDWF which
            # uses PassCalcLocs. For the rest we just specify the previous dir.
            if spawn_count==0:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=False))
            else:
                t.append(CopyVaspOutputs(calc_dir=current_dir, contcar_to_poscar=True))

            t.append(RescaleVolumeTask(initial_pressure=p*1000.0, initial_temperature=1))
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
                                   spawn_count=spawn_count+1,
                                   copy_calcs=copy_calcs,
                                   calc_home=calc_home,
                                   averaging_fraction=averaging_fraction))
            new_fw = Firework(t, name=name)
            return FWAction(stored_data={'pressure': p}, additions=[new_fw])

        else:
            return FWAction(stored_data={'pressure': p, 'density_calculated': True})


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
    required_params = ["calc_home","run_name"]
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
                shutil.copy2(f,target_dir)
            except:
                pass
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
