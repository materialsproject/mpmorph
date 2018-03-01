from pymatgen import Structure, Composition, Element
from pymatgen.io.vasp import Xdatcar
from mpmorph.analysis.clustering_analysis import ClusteringAnalyzer
from mpmorph.analysis.structural_analysis import RadialDistributionFunction
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
import multiprocessing


def process_frame(data):
    frame_num, structure, bond_lengths, prune_els, cluster_els, cluster_bonds = data[0], data[1], data[2], data[3], data[4], data[5]

    # cluster_els = [Element('Si'), Element('O')]
    # cluster_bonds = [('O', 'Si'), ('Si', 'Si')]
    ca = ClusteringAnalyzer(structure, bond_lengths=bond_lengths)
    clusters = ca.get_clusters(prune_els=prune_els, cluster_els=cluster_els, cluster_bonds=cluster_bonds)
    neighbors = ca.cluster_neighbors
    # track_neighbors = ca.track_neighbors
    track_neighbors = []
    del structure, bond_lengths, prune_els
    return frame_num, neighbors, clusters, track_neighbors


class EnvironmentTracker():
    # TODO: Add functionality for multielemental clusters
    def __init__(self):
        pass

    def run(self, structures, cluster_els, cluster_bonds, frames=None, prune_els=[]):
        if frames == None:
            frames = len(structures)
        bond_lens = self.get_bond_distance(structures)

        pool = Pool(multiprocessing.cpu_count(), maxtasksperchild=1000)
        inputs = [(i, structure, bond_lens, prune_els, cluster_els, cluster_bonds) for (i, structure) in enumerate(structures)]
        results = pool.map(process_frame, inputs)
        sort_key = lambda result: result[0]
        sorted(results, key=sort_key)
        neighbor_array = []
        cluster_array = []
        track_neighbors_array = []
        for result in results:
            neighbor_array.append(result[1])
            cluster_array.append(result[2])
            track_neighbors_array.append(result[3])

        pool.close()
        pool.join()
        # neighbor_array = (frames)*[None] #Predeclare 3d array of length of xdatcar
        #         cluster_array = (frames)*[None]
        #         track_neighbor_array = (frames)*[None]
        #         for i in range(frames):
        #             neighbors, clusters, track_neighbors = self.process_frame(structure=structures[len(structures)-frames+i], bond_lengths=bond_lengths, prune_els=prune_els)
        #             track_neighbor_array[i] = track_neighbors
        #             neighbor_array[i] = neighbors
        #             cluster_array[i] = clusters
        return neighbor_array, cluster_array, track_neighbors_array

    def get_bond_distance(self, structures):
        bin_size = 0.01
        cutoff = 5
        rdf = RadialDistributionFunction(structures, step_freq=1, bin_size=bin_size, cutoff=cutoff, smooth=1)
        a = rdf.get_radial_distribution_functions(nproc=multiprocessing.cpu_count())

        coord_nums = []
        bond_lengths = {}
        for pair in a:
            y = a[pair]
            x = np.arange(0, cutoff, bin_size)
            maximum = y[0:int(4 / bin_size)].argmax() * bin_size
            min_past = maximum
            integration_cutoff_i = int(min_past / bin_size) + y[int(min_past / bin_size):].argmin()
            bond_lengths[tuple(sorted(pair))] = integration_cutoff_i * bin_size
        return bond_lengths

    def get_statistics(self, track_el, neighbor_array, structure):
        # How long does it stay near one cluster
        # How many Si does the Li visit?
        # Do Li stick to a cluster or one Si?
        # Do more Si neighbors cause Li to stay longer?
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
                    tracking_list[i][j] += 1

        return tracking_list
