"""
This example generates an AIMD workflow for finding equilibrium density for liquid water at 320 K.
"""

box_scale = 8.9 # edge length of MD box in Angstroms
packmol_path = "~/packmol/packmol/packmol" # Revise as appropriate
# "structure" in this context can be a dict of number of atoms or molecules.
structure = {'H2O':20}
temperature = 320
# Note one can use a pymatgen Structure object also
# E.g. p = Poscar.from_file("POSCAR")
#      structure = p.structure
# MD runs can be backed up in a desired location in the
# order they're run with copy_calcs option:
copy_calcs = True
calc_home = '~/test_H2O_wflows' # Revise as appropriate
# Since we specified a molecule, we must also give the path to xyz
# file of a single sample molecule.
xyz_paths = ['H2O.xyz']


from mpmorph.workflow.workflows import get_wf_density
from fireworks import LaunchPad

amorphous_maker_params = {'box_scale':box_scale, 'packmol_path':packmol_path, 'xyz_paths': xyz_paths}
name = 'H2O_df_'+str(temperature)

wf = get_wf_density(structure, temperature=temperature, nsteps=1000, wall_time=19200, max_rescales=5,
                    amorphous_maker_params=amorphous_maker_params, copy_calcs=copy_calcs, calc_home=calc_home, name=name)

lp = LaunchPad.auto_load()
lp.reset('', require_password=False)
lp.add_wf(wf)