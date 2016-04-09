from pymatgen.io.vasp.sets import MITMDVaspInputSet
from mpmorph.runners.amorphous_maker import AmorphousMaker

# Example:
# LiBO2-V2O5 in a cubic box of edge 11.78 angstroms
a = 11.78
composition = {"V":22, "Li":10, "O":75, "B":10}

glass = AmorphousMaker(composition, a, packmol_path="~/packmol/packmol/packmol")
structure = glass.random_packed_structure


vaspmd = MITMDVaspInputSet(3000, 300, 2000, time_step=2, hubbard_off = False, spin_polarized = True,
           user_incar_settings={'LWAVE':'.False.', 'LCHARG':'.True.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 2, 'ENCUT': 520, 'SIGMA':0.1, 'ALGO':'Very Fast'},
                                settings_file = "MPVaspInputSet.yaml")
vaspmd.write_input(structure, "./new_dir/")