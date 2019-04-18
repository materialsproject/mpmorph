from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen import Structure
from pymatgen.io.vasp import Xdatcar
from itertools import product

import os

__author__ = 'Eric Sivonxay'


class RingsInput:

    @staticmethod
    def write_input(structure, length, filename='input', system_name='vasp_run', time_step=2, user_options={}, bonds={},
                    directory=None):

        if directory is None:
            directory = os.getcwd()

        options ={'real_rdf_discretization':200,
                  'recip_discretization': 500,
                  'max_mod': 25,
                  'smooth_factor': 0.125,
                  'angular_discretization': 90,
                  'real_discretization': 20,
                  'max_depth_rings': 10,
                  'max_depth_chain': 15,
                  'cutoff_radius': 2.6}
        options.update(user_options)

        with open(os.path.join(directory, filename), 'w') as file:
            file.write("#"*39 + "\n")
            file.write("#       R.I.N.G.S. input file         #" + "\n")
            file.write("#"*39 + "\n")
            file.write(system_name + "\n") #1) System Name
            file.write(f'{structure.num_sites} \n') #2) Number of atoms in system
            file.write(f'{len(structure.composition.to_reduced_dict)} \n') #3) Number of chemical species
            elements = " ".join([str(el) for el in structure.composition.elements])
            file.write(elements + "\n") #4) Chemical Species
            file.write(f"{length} \n") #5) Number of MD steps
            file.write("1" + "\n") #6) Lattice type(0: lattice parameters and angles, 1: lattice vectors)
            # 7) Lattice vectors
            for line in structure.lattice.matrix:
                file.write(f'{" ".join([f"{el:<10}" for el in line])} \n')
            # file.write(structure.lattice + "n")
            file.write(f'{time_step} \n') #8) Integration time step of Newton's Equations of Motion
            file.write("VAS" + "\n") #9)  File format for configurations
            file.write("XDATCAR" + "\n") #10) Name of the file
            file.write(f'{options["real_rdf_discretization"]}\n') #11) Real space discretization for the g(r) calculations
            file.write(f'{options["recip_discretization"]}\n') #12) Reciprocal space descretization for the S(q) calculations
            file.write(f'{options["max_mod"]}\n') #13) Maximum modulus of the reciprocal space vectors for the S(q) calculations
            file.write(f'{options["smooth_factor"]}\n') #14) Smoothing factor for the S(q)
            file.write(f'{options["angular_discretization"]}\n') #15) Angular discretization
            file.write(f'{options["real_discretization"]}\n') #16) Real Space discretization for the voids and for the ring stats
            file.write(f'{options["max_depth_rings"]}\n') #17) (Maximum search depth)/2 for ring stats
            file.write(f'{options["max_depth_chain"]}\n') #18) Maximum search depth for chain stats
            file.write('#'*39 + "\n")
            if len(bonds.keys())==0:
                elements = [str(el) for el in structure.composition.elements]
                for i, j in product(elements, elements):
                    pair = f'({str(i)}, {str(j)})'
                    bonds[pair] = 0.0

            for pair in bonds.keys():
                el1 = pair.split(',')[0].replace("(", "")
                el2 = pair.split(',')[1].replace(")", "")
                file.write(f'{el1} {el2}  {bonds[pair]} \n')
            # file.write('{options["angular_discretization"]}\n') #19)
            file.write(f'{options["cutoff_radius"]}\n') #20)
            file.write("#"*39 + "\n")
        return

    @staticmethod
    def write_options(directory=None, user_options={}):
        if directory is None:
            directory = os.getcwd()

        options = {'PBC': True, 'Frac': False, 'g(r)': False,
                   'S(q)': False, 'S(k)': False, 'gfft(r)': False,
                   'MSD': False, 'atMSD': False, 'Bonds': True,
                   'Angles': False}
        chain_options = {'Chains': True, 'Species': 0,
                   'AAAA': False, 'ABAB': False, '1221': False}
        ring_options = {'Rings': True, 'Species': 0, 'ABAB': False,
                        'Rings0': False, 'Rings1': False, 'Rings2': False,
                        'Rings3': True, 'Rings4': True, 'Prim_Rings': True,
                        'Str_Rings': False, 'BarycRings': False,
                        'Prop-1': False, 'Prop-2': False, 'Prop-3': False,
                        'Prop-4': False, 'Prop-5': False}
        other_options = {"Vacuum": False, 'Evol': False, 'Dxout': False,
                         "RadOut": False, "RingsOut": False, "DRngOut": False,
                         "VoidsOut": False, 'TetraOut': False, 'TrajOut': False,
                         "Output": "my-output.out"}
        options.update(user_options)
        with open(os.path.join(directory, 'options'), 'w') as file:
            file.write("#" * 39 + "\n")
            file.write("#       R.I.N.G.S. options file         #" + "\n")
            file.write("#" * 39 + "\n")
            for key, item in options.items():
                if type(item) is bool:
                    val = '.true.' if item else '.false.'
                    file.write(f'{key:<20}{val}\n')
                else:
                    file.write(f'{key:<20}{item:^7}\n')
            file.write('-'*39 + '\n')
            for key, item in chain_options.items():
                if type(item) is bool:
                    val = '.true.' if item else '.false.'
                    file.write(f'{key:<20}{val}\n')
                else:
                    file.write(f'{key:<20}{item:^7}\n')

            file.write('-'*39 + '\n')
            for key, item in ring_options.items():
                if type(item) is bool:
                    val = '.true.' if item else '.false.'
                    file.write(f'{key:<20}{val}\n')
                else:
                    file.write(f'{key:<20}{item:^7}\n')

            file.write('-'*39 + '\n')
            for key, item in other_options.items():
                if type(item) is bool:
                    val = '.true.' if item else '.false.'
                    file.write(f'{key:<20}{val}\n')
                else:
                    file.write(f'{key:<20}{item:^7}\n')
            file.write("#" * 39)
        return


# TODO: add ability to handle multiple cutoff radii (one for each pair)
# class Rings(object):
#     #
#     def __init__(self, directory, max_ring_size = 5, radius=2.6):
#         self.directory = directory or os.getpwd()
#         self.radius = radius
#
#     def process_md(self, xdatcar):
#         structures = xdatcar.structures
#         snap_rings = []
#         #for each snapshot
#         #for _structure in structures:
#         _structure = structures[0]
#         snap_rings.append(self.find_rings(_structure))
#         return snap_rings
#
#     def find_rings(self, structure):
#         struct_sites = structure.sites
#
#         distance_matrix = self.get_distance_matrix(structure)
#         neighbors_matrix = self.get_neighbors(distance_matrix=distance_matrix, structure=structure)
#         # start at each site
#         # Potential to parallelize
#
#         for i in range(len(struct_sites)):
#             ring_size = 0
#             _neighbors = neighbors_matrix[i]
#
#             if len() > max_ring_size:
#                 break
#         return
#
#     def ring_search(self, neighbors_matrix, struct_sites, start_site, neighbor_sites, current_site, ring_nodes):
#
#         self.ring_search(neighbors_matrix=neighbors_matrix, struct_sites=struct_sites, current_site=current_site, ring_nodes=ring_nodes)
#         return
#
#     def get_distance_matrix(self, structure):
#         struct_sites = structure.sites
#         distance_matrix = [[0 for x in struct_sites] for y in struct_sites]
#         for i in range(len(struct_sites)):
#             for j in range(i):
#                 _distance = struct_sites[i].distance(struct_sites[j])
#                 distance_matrix[j][i] = _distance
#                 distance_matrix[i][j] = _distance
#         return distance_matrix
#
#     def get_neighbors(self, distance_matrix, structure):
#         neighbors = [[] for x in self.struct_sites]
#         for i in range(len(self.struct_sites)):
#             for j in range(len(self.struct_sites)):
#                 if distance_matrix[i][j] <= self.radius:
#                     neighbors[i].append(j)
#         return neighbors
#
#     def import_data(self):
#         pass
#
#     def write_options(self):
#         pass
#
#     def write_input(self, xdatcar, filename='input', system_name='vasp_run', time_step=2, user_options={}, bonds={}):
#         options ={'real_rdf_discretization':200,
#                   'recip_discretization': 3,
#                   'max_mod': 3,
#                   'smooth_factor': 1,
#                   'angular_discretization': 13,
#                   'real_discretization': 10,
#                   'max_depth_rings': 2,
#                   'max_depth_chain': 2,
#                   'cutoff_radius': 10}
#         options.update(user_options)
#         structures = xdatcar.structures
#         _structure = structures[0]
#
#         with open(os.path.join(self.directory, filename), 'w') as file:
#             file.write("#"*39 + "\n")
#             file.write("#       R.I.N.G.S. input file         #" + "\n")
#             file.write("#"*39 + "\n")
#             file.write(system_name + "\n") #1) System Name
#             file.write(_structure.num_sites) #2) Number of atoms in system
#             file.write(len(_structure.composition.to_reduced_dict)) #3) Number of chemical species
#             elements = " ".join([str(el) for el in _structure.composition.elements])
#             file.write(elements + "\n") #4) Chemical Species
#             file.write(len(structures) + "\n") #5) Number of MD steps
#             file.write("1" + "\n") #6) Lattice type(0: lattice parameters and angles, 1: lattice vectors)
#             # 7) Lattice vectors
#             for line in _structure.lattice.matrix:
#                 file.write(f'{" ".join(line)} \n')
#             file.write(_structure.lattice + "n")
#             file.write(time_step + "\n") #8) Integration time step of Newton's Equations of Motion
#             file.write("VAS" + "\n") #9)  File format for configurations
#             file.write("XDATCAR" + "\n") #10) Name of the file
#             file.write('{options["real_rdf_discretization"]}\n') #11) Real space discretization for the g(r) calculations
#             file.write('{options["recip_discretization"]}\n') #12) Reciprocal space descretization for the S(q) calculations
#             file.write('{options["max_mod"]}\n') #13) Maximum modulus of the reciprocal space vectors for the S(q) calculations
#             file.write('{options["smooth_factor"]}\n') #14) Smoothing factor for the S(q)
#             file.write('{options["angular_discretization"]}\n') #15) Angular discretization
#             file.write('{options["real_discretization"]}\n') #16) Real Space discretization for the voids and for the ring stats
#             file.write('{options["max_depth_rings"]}\n') #17) (Maximum search depth)/2 for ring stats
#             file.write('{options["max_depth_chain"]}\n') #18) Maximum search depth for chain stats
#             file.write('#'*39 + "\n")
#             for pair in bonds.keys():
#                 file.write(f'{pair[0]} {pair[1]}  {bonds[pair]} \n')
#             # file.write('{options["angular_discretization"]}\n') #19)
#             file.write('{options["cutoff_radius"]}\n') #20)
#             file.write("#"*39 + "\n")
#         return
#
#     def run_rings(self):
#         pass