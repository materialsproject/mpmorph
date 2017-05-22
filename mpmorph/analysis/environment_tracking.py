from pymatgen import Structure, Composition, Element
from pymatgen.io.vasp import Xdatcar
from mpmorph.analysis.clustering_analysis import ClusteringAnalyzer
from mpmorph.analysis.structural_analysis import RadialDistributionFunction
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

class EnvironmentTracker():
    # TODO: Add functionality for multielemental clusters
    def __init__(self):
        pass

    def run(self, structures, frames=None, prune_els=[]):
        if frames==None:
            frames = len(structures)
        bond_lengths = self.get_bond_distance(structures)

        neighbor_array = (frames)*[None] #Predeclare 3d array of length of xdatcar
        cluster_array = (frames)*[None]
        track_neighbor_array = (frames)*[None]
        for i in range(frames):
            neighbors, clusters, track_neighbors = self.process_frame(structure=structures[len(structures)-frames+i], bond_lengths=bond_lengths, prune_els=prune_els)
            track_neighbor_array[i] = track_neighbors
            neighbor_array[i] = neighbors
            cluster_array[i] = clusters
        return neighbor_array, cluster_array, track_neighbor_array

    def process_frame(self, structure, bond_lengths, prune_els=[]):
        ca = ClusteringAnalyzer(structure, bond_lengths=bond_lengths)
        clusters = ca.get_clusters(prune_els=prune_els)
        neighbors = ca.cluster_neighbors
        track_neighbors = ca.track_neighbors
        return neighbors, clusters, track_neighbors

    def get_bond_distance(self, structures):
        bin_size = 0.1
        cutoff = 5
        rdf = RadialDistributionFunction(structures, step_freq=1, bin_size=bin_size, cutoff=cutoff, smooth=1)
        a = rdf.get_radial_distribution_functions(nproc=multiprocessing.cpu_count())
        rdf.plot_radial_distribution_functions()
        plt.show()

        coord_nums = []
        bond_lengths = {}
        for pair in a:
            y = a[pair]
            x = np.arange(0, cutoff, bin_size)
            maximum = y[0:int(4 / bin_size)].argmax() * bin_size
            min_past = maximum
            integration_cutoff_i = int(min_past / bin_size) + y[int(min_past / bin_size):].argmin()
            bond_lengths[pair] = integration_cutoff_i*bin_size
        return bond_lengths

    def get_statistics(self, track_el, neighbor_array, structure):
        #How long does it stay near one cluster
        #How many Si does the Li visit?
        #Do Li stick to a cluster or one Si?
        #Do more Si neighbors cause Li to stay longer?
        struct_sites = structure.sites
        track_positions = []
        cluster_positions = []
        for i in range(len(structure.species)):
            if structure.species[i] == track_el:
                track_positions.append(i)
            else:
                cluster_positions.append(i)
        tracking_list = [[0 for x in range(len(cluster_positions))] for y in range(len(track_positions))]

        for frame in neighbor_array:
            for i in range(len(track_positions)):
                for j in frame[i]:
                    tracking_list[i][j]+=1

        return tracking_list
