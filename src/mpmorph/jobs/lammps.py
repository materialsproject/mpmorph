from subprocess import PIPE, Popen
from jobflow import Maker, job
from pymatgen.io.lammps.utils import LammpsRunner
import logging

from typing import List

from pymatgen.io.lammps.inputs import LammpsTemplateGen
from pymatgen.io.lammps.outputs import LammpsDump, parse_lammps_log
from pymatgen.io.lammps.data import LammpsData
from pymatgen.core.structure import Structure

class RunLammpsMaker(Maker):
    """
    Run LAMMPS directly (no custodian).
    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    name = "RUN_LAMMPS"

    @job
    def make(self, lammps_bin: str,
                   script_template_path: str,
                   script_options: dict,
                   structure: Structure = None,
                   log_filename: str = "log.lammps",
                   data_filename: str = "data.lammps",
                   dump_files: List[str] = None):

        data = LammpsData.from_structure(structure, atom_style='atomic')
        # Write the input files
        linp = LammpsTemplateGen().get_input_set(script_template=script_template_path,
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

        dump_files = dump_files or []
        dump_files = [dump_files] if isinstance(dump_files, str) else dump_files

        # Construct various dumps objects
        dumps = []
        if dump_files:
            for df in dump_files:
                dumps.append((df, LammpsDump.from_file(df)))

        # Construct log object
        log = parse_lammps_log(log_filename)

        return {
            "dumps": dumps,
            "log": log
        }