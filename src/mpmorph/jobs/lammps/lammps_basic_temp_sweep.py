import os
from jobflow import Maker, job

import pandas as pd
from pymatgen.core.structure import Structure

from mpmorph.schemas.lammps_calc import LammpsCalc
from .helpers import run_lammps, trajectory_from_lammps_dump

from pkg_resources import resource_filename

class BasicLammpsTempSweepMaker(Maker):
    """
    Run LAMMPS directly using m3gnet sweeping over a range of temperatures.
    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    name = "LAMMPS_CALCULATION"

    @job(trajectory="trajectory", output_schema=LammpsCalc)
    def make(self, temp_initial: int,
                   temp_final: int,
                   total_steps: int,
                   structure: Structure = None):

        lammps_bin = os.environ.get("LAMMPS_CMD")
        m3gnet_path = os.environ.get("M3GNET_PATH")

        chem_sys_str = " ".join(el.symbol for el in structure.composition.elements)

        script_options = {
            "tempstart": temp_initial,
            "tempstop": temp_final,
            "m3gnet_path": m3gnet_path,
            "species": chem_sys_str,
            "total_steps": total_steps,
            "print_every_n_step": 10
        }
        

        template_path = resource_filename('mpmorph', 'jobs/lammps/templates/basic_temp_sweep.lammps')

        run_lammps(structure, template_path, script_options, lammps_bin)

        trajectory = trajectory_from_lammps_dump("trajectory.lammpstrj")

        # df = pd.read_csv("step_temp_vol_density.txt", delimiter=" ", index_col="step", skiprows=1, names=["step", "temp", "vol", "density"])

        metadata = {
            "temp_initial": temp_initial,
            "temp_final": temp_final,
            "total_steps": total_steps
        }

        output = LammpsCalc(
            dir_name=os.getcwd(),
            trajectory=trajectory,
            composition=structure.composition,
            reduced_formula=structure.composition.reduced_formula,
            metadata=metadata,
            dump_data={}
        )
        return output