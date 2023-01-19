from jobflow import Maker, job
from pymatgen.io.lammps.utils import LammpsRunner
from pymatgen.io.lammps.inputs import La
import logging
import os

from pymatgen.io.lammps.inputs import TemplateInputGen, LammpsTemplateGen
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

    @job
    def make(self, lammps_cmd: str,
                   script_template_path: str,
                   script_options: dict,
                   structure: Structure = None,
                   data_filename: str = None,
                   log_filename: str = "log.lammps"):

        if data_filename is None:
            data = LammpsData.from_structure(structure)
        else:
            data = LammpsData.from_file(data_filename)
        # Write the input files
        linp = LammpsTemplateGen().get_input_set(script_template=script_template_path,
                                                 settings=script_options, 
                                                 data=data,
                                                 data_filename='data.dump')

        input_name = f'lammps.in'
        linp.write_input(input_name)

        # Run LAMMPS
        lmps_runner = LammpsRunner(input_name, lammps_cmd)
        stdout, stderr = lmps_runner.run()
        logging.info(f"LAMMPS finished running: {stdout} \n {stderr}")

        dump_files = dump_files or []
        dump_files = [dump_files] if isinstance(dump_files, str) else dump_files

        # Construct various dumps objects
        dumps = []
        if dump_files:
            for df in dump_files:
                dumps.append((df, LammpsDump.from_file(df)))

        # Construct log object
        log = parse_lammps_log(log_filename)

        logging.info(f"Getting task doc for base dir :{path}")
        d = self.generate_doc(path, lmps_input, log, dumps)        