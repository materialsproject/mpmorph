from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen import Structure
from pymatgen.io.vasp import Xdatcar

import os

class Rings(object):
    def __init__(self, directory):
        self.directory = directory or os.getpwd()

    def import_data(self):
        pass

    def write_input(self):
        pass

    def write_options(self, xdatcar, filename='input', system_name='vasp_run', time_step=2):
        structures = xdatcar.structures
        _structure = structures[0]

        file = open(os.path.join(self.directory, filename), 'w')
        file.write("#"*39 + "\n")
        file.write("#       R.I.N.G.S. input file         #" + "\n")
        file.write("#"*39 + "\n")
        file.write(system_name + "\n") #1) System Name
        file.write(_structure.num_sites) #2) Number of atoms in system
        file.write(len(_structure.composition.to_reduced_dict)) #3) Number of chemical species
        elements = ""
        for el in _structure.composition.elements:
            elements = elements + str(el) + " "
        file.write(elements + "\n") #4) Chemical Species
        file.write(len(structures) + "\n") #5) Number of MD steps
        file.write("1" + "\n") #6) Lattice type(0: lattice parameters and angles, 1: lattice vectors)
        file.write(_structure.lattice + "n") #7) Lattice vectors
        file.write(time_step + "\n") #8) Integration time step of Newton's Equations of Motion
        file.write("VAS" + "\n") #9)  File format for configurations
        file.write("XDATCAR" + "\n") #10) Name of the file
        file.write( + "\n") #11) Real space discretization for the g(r) calculations
        file.write( + "\n") #12) Reciprocal space descretization for the S(q) calculations
        file.write( + "\n") #13) Maximum modulus of the reciprocal space vectors for the S(q) calculations
        file.write( + "\n") #14) Smoothing factor for the S(q)
        file.write( + "\n") #15)
        file.write( + "\n") #16) Real Space discretization for the voids and for the ring stats
        file.write( + "\n") #17) (Maximum search depth)/2 for ring stats
        file.write( + "\n") #18) Maximum search depth for chain stats
        file.write( + "\n") #19)
        file.write( + "\n") #20)
        file.write("#"*39 + "\n")
        file.write("#"*39 + "\n")
        return

    def run_rings(self):
        pass