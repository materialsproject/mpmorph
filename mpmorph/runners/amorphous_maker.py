from __future__ import division
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from collections import OrderedDict
from pymatgen.io.vasp.sets import MITMDVaspInputSet
import numpy as np
import os


class AmorphousMaker(object):
    def __init__(self, el_num_dict, box_scale, tol = 2.0, packmol_path="packmol"):
        """
        Class for generating initial constrained-random packed structures for the
        simulation of amorphous or liquid structures. This is a wrapper for "packmol".
        Only works for cubic boxes for now.

        Args:
        el_num_dict: dictionary of number of atoms of each species
            e.g. {"V":22, "Li":10, "O":75, "B":10}
        box_scale: all lattice vectors are multiplied with this scalar value
            e.g. edge length of a cubic simulation box
        """
        self.el_num_dict = el_num_dict
        self.box_scale = box_scale
        self.tol = tol
        self._lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self._structure = None
        self._el_dict = None
        self.packmol_path = packmol_path

    def __repr__(self):
        return "AmorphousMaker: generates constrained-random packed initial structure for MD."

    @property
    def box(self):
        """
        Returns: box vectors scaled with box_scale
        """
        return (np.array(self._lattice)*self.box_scale).tolist()

    @property
    def random_packed_structure(self):
        """
        Returns:
        """
        self._el_dict = self.call_packmol(clean=True)
        self._structure = self.get_structure(self._el_dict, self.box)
        return self._structure

    def call_packmol(self, clean):
        """
        Args:
            clean: (bool) remove intermediate files during generation
        Returns:
            dict of coordinates
        """

        # this ensures periodic boundaries don't cause problems
        pm_l = self.tol/2
        pm_h = self.box_scale-self.tol/2

        with open("packmol.input", "w") as f:
            f.write("tolerance "+ str(self.tol) +"\nfiletype xyz\noutput mixture.xyz\n")
            for el in self.el_num_dict:
                f.write("structure " + el + ".xyz\n" + "  number " + str(self.el_num_dict[el])
                          + "\n  inside box" + 3*(" " + str(pm_l)) + 3*(" " + str(pm_h))
                          + "\nend structure\n\n")
        for el in self.el_num_dict.keys():
           with open(el+".xyz", "w") as f:
               f.write("1\ncomment\n" + el +" 0.0 0.0 0.0\n")
        try:
            os.system(self.packmol_path + " < packmol.input")
        except:
            raise OSError("packmol cannot be found!")
        if clean:
            for el in el_num_dict.keys():
                os.system("rm " + el + ".xyz")
            os.system("rm packmol.input")
        return self.xyz2dict("mixture.xyz", clean=clean)

    def xyz2dict(self, filename, clean=False):
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
            raise ValueError("Inconsistent number of atoms")
        self._el_dict=OrderedDict(el_dict)
        if clean:
            os.system("rm " + filename)
        return self._el_dict

    @staticmethod
    def get_structure(el_dict, lattice):
        """
        Args:
            el_dict: is a dictionary of coordinates of each species
                e.g. {'Cr': [ [0.0, 0.0, 0.0], [1.2, 1.2, 1.2], ...], 'O': [...] }
            lattice: is the lattice in the form of [[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]
        Returns: pymatgen Structure
        """
        species = []
        coords = []
        for el in el_dict.keys():
            for atom in el_dict[el]:
                species.append(el)
                coords.append(atom)
        return Structure(lattice, species, coords, coords_are_cartesian=True)

    def get_poscar(self):
        return Poscar(self.random_packed_structure)


    @staticmethod
    def xyzdict_to_poscar(el_dict, lattice, filepath ="POSCAR"):
        """
        Generates XYZ file from element coordinate dictionary and lattice
        Args:
            el_dict: is a dictionary of coordinates of each species
            lattice: is the lattice in the form of [[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]
            filepath: path to POSCAR to be generated
        Returns:
            writes a POSCAR file
        """
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


if __name__ == "__main__":

    # Example:
    # LiBO2-V2O5 in a cubic box of edge 11.78 angstroms
    a = 11.78
    el_num_dict = {"V":22, "Li":10, "O":75, "B":10}

    glass = AmorphousMaker(el_num_dict, a, packmol_path="~/packmol/packmol/packmol")
    struct = glass.random_packed_structure


    vaspmd = MITMDVaspInputSet(3000, 300, 2000, time_step=2, hubbard_off = False, spin_polarized = True,
           user_incar_settings={'LWAVE':'.False.', 'LCHARG':'.True.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 2, 'ENCUT': 520, 'SIGMA':0.1, 'ALGO':'Very Fast'},
                                settings_file = "MPVaspInputSet.yaml")
    vaspmd.write_input(struct, "./new_dir/")





