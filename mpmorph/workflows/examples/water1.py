"""
This example generates a dynamic AIMD workflow for finding equilibrium density for liquid water at 320 K.
"""

# Parameters:
box_scale = 8.9  # edge length of MD box in Angstroms, can also be a numpy array that scales the lattice
packmol_path = "~/packmol/packmol/packmol"  # Revise as appropriate
structure = {'H2O': 20}  # "structure" in this context can be a dict of number of atoms or molecules.
temperature = 320

# Note one can use a pymatgen Structure object also
# E.g. p = Poscar.from_file("POSCAR")
#      structure = p.structure

copy_calcs = True  # MD runs can be backed up in a desired location
calc_home = '~/test_H2O_wflows'  # This is the location to copy the calculations if copy_calcs=True

# Since we specified a molecule, we must also give the path to xyz
# file of a single sample molecule.
xyz_paths = ['H2O.xyz']
name = 'H2O_df_' + str(temperature)

from fireworks import LaunchPad
from mpmorph.workflows.old_workflows import get_wf_density

amorphous_maker_params = {'box_scale': box_scale, 'packmol_path': packmol_path, 'xyz_paths': xyz_paths, 'tol': 2.0}

wf = get_wf_density(structure, temperature=temperature, pressure_threshold=0.5, nsteps=1000, wall_time=19200,
                    max_rescales=5,
                    amorphous_maker_params=amorphous_maker_params, copy_calcs=copy_calcs, calc_home=calc_home,
                    name=name)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
