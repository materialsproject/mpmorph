from subprocess import PIPE, Popen

import os
from jobflow import Maker, job

import pandas as pd

from pymatgen.io.lammps.inputs import LammpsTemplateGen
from pymatgen.io.lammps.data import LammpsData
from pymatgen.core.structure import Structure

from mpmorph.schemas.pv_data_doc import MDPVDataDoc

from pkg_resources import resource_filename

class LammpsVolMaker(Maker):
    """
    Run LAMMPS directly using m3gnet (no custodian).
    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    name = "LAMMPS_TO_VOLUME"

    @job
    def make(self, temperature: int,
                   total_steps: int,
                   structure: Structure = None):

        lammps_bin = os.environ.get("LAMMPS_CMD")
        m3gnet_path = os.environ.get("M3GNET_PATH")

        chem_sys_str = " ".join(el.symbol for el in structure.composition.elements)
        script_options = {
            "temperature": temperature,
            "m3gnet_path": m3gnet_path,
            "species": chem_sys_str,
            "total_steps": total_steps,
            "print_every_n_step": 10
        }                   
        

        template_path = resource_filename('mpmorph', 'jobs/lammps-templates/template.lammps')

        data_filename: str = "data.lammps"
        data = LammpsData.from_structure(structure, atom_style='atomic')
        # Write the input files
        linp = LammpsTemplateGen().get_input_set(script_template=template_path,
                                                 settings=script_options, 
                                                 data=data,
                                                 data_filename=data_filename)

        linp.write_input(directory=".")
        input_name = "in.lammps"
        # Run LAMMPS

        lammps_cmd = [lammps_bin, "-in", input_name]
        print(f"Running: {' '.join(lammps_cmd)}")
        with Popen(lammps_cmd, stdout=PIPE, stderr=PIPE) as p:
            (stdout, stderr) = p.communicate()

        print(f"LAMMPS finished running: {stdout} \n {stderr}")


        avging_window = int(min(total_steps / 100, 50))
        df = pd.read_csv("step_temp_vol_density.txt", delimiter=" ", skiprows=1, names=["step", "temp", "vol", "density"])
        eq_vol = df.iloc[-avging_window::]['vol'].values.mean()

        return eq_vol