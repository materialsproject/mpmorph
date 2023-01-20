from subprocess import PIPE, Popen
from jobflow import Maker, job

import pandas as pd

from pymatgen.io.lammps.inputs import LammpsTemplateGen
from pymatgen.io.lammps.data import LammpsData
from pymatgen.core.structure import Structure

from mpmorph.schemas.pv_data_doc import MDPVDataDoc

from pkg_resources import resource_filename

class LammpsVolMaker(Maker):
    """
    Run LAMMPS directly (no custodian).
    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    name = "RUN_LAMMPS"

    @job
    def make(self, lammps_bin: str,
                   script_options: dict,
                   structure: Structure = None,
                   data_filename: str = "data.lammps"):
        
        template_path = resource_filename('mpmorph', 'jobs/lammps-templates/template.lammps')


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


        filecontents = pd.read_csv("step_temp_vol_density.txt", delimiter=" ", skiprows=1, index_col="step", names=["step", "temp", "vol", "density"])
        eq_vol = filecontents[["vol"]].iloc[-1].values[0]

        return eq_vol