from ase.io.lammpsrun import read_lammps_dump_text
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.lammps.inputs import LammpsTemplateGen
from pymatgen.io.lammps.data import LammpsData
from subprocess import PIPE, Popen

def trajectory_from_lammps_dump(dump_path):
    with open(dump_path, "r+") as f:
        atoms = read_lammps_dump_text(f, index=slice(0, None))

    structs = []

    for a in atoms:
        structs.append(AseAtomsAdaptor().get_structure(a))

    return Trajectory.from_structures(structs, constant_lattice=False)

def run_lammps(structure, template_path, template_opts, lammps_bin):
    data_filename: str = "data.lammps"
    data = LammpsData.from_structure(structure, atom_style='atomic')
    # Write the input files
    linp = LammpsTemplateGen().get_input_set(script_template=template_path,
                                                settings=template_opts, 
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