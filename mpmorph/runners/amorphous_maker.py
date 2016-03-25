from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from collections import OrderedDict
from pymatgen.io.vasp.sets import MITMDVaspInputSet
import sys
import os

def xyz2dict(filename):
    """
    This is a generic xyz to dictionary convertor. 
    Used to get the structure from packmol output.
    """ 
    with open(filename,'r') as f:
        lines = f.readlines()
        N = int(lines[0].rstrip('\n'))
        el_dict = {}
        for line in lines[2:]:
            l = line.rstrip('\n').split()
            if l[0] in el_dict:
                el_dict[ l[0] ].append( [ float (i) for i in l[1:] ] )
            else:
                el_dict[ l[0] ] = [ [ float (i) for i in l[1:] ] ] 
    if N != sum( [ len(x) for x in el_dict.values() ] ):
        raise "Inconsistent number of atoms" 
    return OrderedDict(el_dict)

def get_structure(el_dict, lattice):
    """
    Args: el_dict is a dictionary of coordinates of each species
          lattice is the lattice in the form of [[],[],[]]
    Returns: pymatgen Structure
    """
    species = []
    coords = []
    for el in el_dict.keys():
        for atom in el_dict[el]:
            species.append(el)
            coords.append(atom)
    struct = Structure(lattice, species, coords, coords_are_cartesian = True)
    return struct 

def xyzdict_to_poscar(el_dict, lattice, filepath ="POSCAR"):
    with open(filepath, "w") as f:
        f.write("Parsed form XYZ file\n")
        f.write("1.0\n")
        for vec in lattice:
            f.write(" ".join( [str(v) for v in vec ] ) + "\n")
        el_dict = OrderedDict(el_dict)
        for key in el_dict.keys():
            f.write(key + " ")
        f.write("\n")
        for key in el_dict.keys():
            f.write( str(len(el_dict[key])) + " ")
        f.write("\nCartesian\n")
        for key in el_dict.keys():
            for atom in el_dict[key]:
                f.write( " ".join( [ str(i) for i in atom ]) + "\n")

def call_packmol(el_num_dict, box_size, tol = 2.0, packmol_path = "packmol"):
    """
    Args:
        - el_num_dict: dictionary of number of atoms of each species
        - tol: hard-sphere radius in angstorms
        - lattice: lattice vectors in the form of [[],[],[]]
    Returns:
        packmol writes an xyz file named mixture.xyz
    """
    # this ensures periodic boundaries don't cause problems
    pm_l = tol/2
    pm_h = box_size-tol/2

    with open("packmol.input", "w") as f:
        f.write("tolerance "+ str(tol) +"\nfiletype xyz\noutput mixture.xyz\n")
        for el in el_num_dict:
            f.write("structure " + el + ".xyz\n" + "  number " + str(el_num_dict[el])
                      + "\n  inside box" + 3*(" " + str(pm_l)) + 3*(" " + str(pm_h))
                      + "\nend structure\n\n") 
    for el in el_num_dict.keys():
       with open(el+".xyz", "w") as f:
           f.write("1\ncomment\n" + el +" 0.0 0.0 0.0\n")
    try:
        os.system(packmol_path + " < packmol.input")
    except:
        raise "packmol could not be found!" 
    for el in el_num_dict.keys():
        os.system("rm " + el + ".xyz")
    os.system("rm packmol.input")
    return xyz2dict("mixture.xyz")
   
if __name__ == "__main__":
     # Testing functions:
     a = 11.78
     el_num_dict = {"V":22, "Li":10, "O":75, "B":10}     

     el_dict = call_packmol(el_num_dict, 11.78, 
               tol = 2.0, packmol_path = "~/packmol/packmol/packmol")

     lattice = [ [a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a] ]
     xyzdict_to_poscar(el_dict, lattice, filepath = "POSCAR_xyzdict_toposcar")
     struct = get_structure(el_dict, lattice)
     poscar = Poscar(struct)
     poscar.write_file("POSCAR_from_struct")
     
     vaspmd = MITMDVaspInputSet(3000, 300, 2000, time_step=2, hubbard_off = False, spin_polarized = True,
           user_incar_settings={'LWAVE':'.False.', 'LCHARG':'.True.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 2, 'ENCUT': 520, 'SIGMA':0.1, 'ALGO':'Very Fast'},
                                settings_file = "MPVaspInputSet.yaml")
     vaspmd.write_input(struct, "./new_dir/")


     
     
    
