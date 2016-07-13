from __future__ import division
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from collections import OrderedDict
import numpy as np
import os
import shutil


class AmorphousMaker(object):
    def __init__(self, el_num_dict, box_scale, tol=2.0, packmol_path="packmol", clean=True, xyz_paths=None):
        """
        Class for generating initial constrained-random packed structures for the
        simulation of amorphous or liquid structures. This is a wrapper for "packmol" package.
        Only works for cubic boxes for now.
        Args:
            el_num_dict (dict): dictionary of number of atoms of each species. If
                number of molecules is specified, an xyz file with the same name needs to be provided as xyz_paths.
                e.g. {"V":22, "Li":10, "O":75, "B":10}
                e.g. {"H2O": 20}
            box_scale (float): all lattice vectors are multiplied with this scalar value.
                e.g. edge length of a cubic simulation box
            tol (float): tolerance factor for how close the atoms can get (angstroms).
                e.g. tol = 2.0 angstroms
            packmol_path (str): path to the packmol executable
            clean (bool): whether the intermedite files generated are deleted.
            xyz_paths (list): list of paths (str) to xyz files correpsonding to molecules, if given so in el_num_dict.
                file names must match the molecule formula.
        """
        self.el_num_dict = el_num_dict
        self.box_scale = box_scale
        self.tol = tol
        self._lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self._structure = None
        self._el_dict = None
        self.packmol_path = packmol_path
        self.clean = clean
        self.xyz_paths = xyz_paths
        if self.xyz_paths:
            assert len(self.xyz_paths)==len(self.el_num_dict.keys())
            self.clean = False

    def __repr__(self):
        return "AmorphousMaker: generates constrained-random packed initial structure for MD using packmol."

    @property
    def box(self):
        """
        Returns: box vectors scaled with box_scale
        """
        return (np.array(self._lattice)*self.box_scale).tolist()

    @property
    def random_packed_structure(self):
        """
        Returns: A constrained-random packed Structure object
        """
        self._el_dict = self.call_packmol()
        self._structure = self.get_structure(self._el_dict, self.box)
        return self._structure

    def call_packmol(self):
        """
        Returns:
            A dict of coordinates of atoms for each element type
            e.g. {'V': [[4.969925, 8.409291, 5.462153], [9.338829, 9.638388, 9.179811], ...]
                  'Li': [[5.244308, 8.918049, 1.014577], [2.832759, 3.605796, 2.330589], ...]}
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

            if self.xyz_paths:
                for path in self.xyz_paths:
                    try:
                        shutil.copy2(path,'./')
                    except:
                        pass
            else:
                for el in self.el_num_dict.keys():
                    with open(el+".xyz", "w") as f:
                       f.write("1\ncomment\n" + el + " 0.0 0.0 0.0\n")
        try:
            os.system(self.packmol_path + " < packmol.input")
        except:
            raise OSError("packmol cannot be found!")
        if self.clean:
            for el in self.el_num_dict.keys():
                os.system("rm " + el + ".xyz")
            os.system("rm packmol.input")
        return self.xyz_to_dict("mixture.xyz")

    def xyz_to_dict(self, filename):
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
        if self.clean:
            os.system("rm " + filename)
        return self._el_dict

    @staticmethod
    def get_structure(el_dict, lattice):
        """
        Args:
            el_dict (dict): coordinates of atoms for each element type
            e.g. {'V': [[4.969925, 8.409291, 5.462153], [9.338829, 9.638388, 9.179811], ...]
                  'Li': [[5.244308, 8.918049, 1.014577], [2.832759, 3.605796, 2.330589], ...]}
            lattice (list): is the lattice in the form of [[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]
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
            el_dict (dict): coordinates of atoms for each element type
            e.g. {'V': [[4.969925, 8.409291, 5.462153], [9.338829, 9.638388, 9.179811], ...]
                  'Li': [[5.244308, 8.918049, 1.014577], [2.832759, 3.605796, 2.330589], ...]}
            lattice (list): is the lattice in the form of [[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]
            filepath (str): path to POSCAR to be generated
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