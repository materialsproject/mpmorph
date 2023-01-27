from subprocess import PIPE, Popen

import os
from jobflow import Maker, job

import pandas as pd

from pymatgen.io.lammps.inputs import LammpsTemplateGen
from pymatgen.io.lammps.data import LammpsData
from pymatgen.core.structure import Structure
from ase.io.lammpsrun import read_lammps_dump_text
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.trajectory import Trajectory

from mpmorph.schemas.lammps_calc import LammpsCalc

from pkg_resources import resource_filename

class LammpsCalcMaker(Maker):
    """
    Run LAMMPS directly using m3gnet (no custodian).
    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    name = "LAMMPS_CALCULATION"

    @job(trajectory="trajectory", output_schema=LammpsCalc)
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

        # Build trajectory from LAMMPS output .xyz file
        with open("trajectory.lammpstrj", "r+") as f:
            atoms = read_lammps_dump_text(f, index=slice(0, None))

        structs = []

        for a in atoms:
            structs.append(AseAtomsAdaptor().get_structure(a))

        trajectory = Trajectory.from_structures(structs, constant_lattice=False)

        df = pd.read_csv("step_temp_vol_density.txt", delimiter=" ", index_col="step", skiprows=1, names=["step", "temp", "vol", "density"])

        metadata = {
            "temperature": temperature,
            "total_steps": total_steps
        }

        output = LammpsCalc(
            dir_name=os.getcwd(),
            trajectory=trajectory,
            composition=structure.composition,
            reduced_formula=structure.composition.reduced_formula,
            metadata=metadata,
            dump_data=df.to_dict()
        )
        return output