from pymatgen import Site, Structure, Element
import numpy as np

class ClusteringAnalyzer(object):
    """

    """

    def __init__(self, input_structure, bond_lengths):
        self.bond_lengths = bond_lengths
        self.radius = 1000
        self.input_structure = input_structure
        self.clusters = []
        self.struct_sites = []
        self.track_distance_matrix = []
        self.cluster_distance_matrix = []
        self.neighbors = []

    def get_clusters(self, cluster_els=[Element("Si")], track_els=[Element("Li")], prune_els = None):
        '''
        TODO: Add functionality to remove elements other than the desired
              RDF's of Cluster
              Return Sites in Cluster
        :param prune_els: List of Elements to exclude from cluster analysis
        :return: clusters
        '''
        #Remove all elements that are not tracked or part of the clusters
        if prune_els != None:
            pruned_structure = self.input_structure.copy()
            pruned_structure.remove_species(prune_els)
        self.radius = self.bond_lengths[(cluster_els[0].name, cluster_els[0].name)]

        #Get structure with only cluster elements
        cluster_structure = self.input_structure.copy()
        els_in_structure = cluster_structure.composition.elements
        for el in els_in_structure:
            if el not in cluster_els:
                cluster_structure.remove_species([el])

        self.track_distance_matrix = self.get_track_distance_matrix(structure=pruned_structure, track_els=track_els)
        self.cluster_distance_matrix = self.get_distance_matrix(structure=cluster_structure)
        self.track_neighbors = self.get_neighbors(distance_matrix = self.track_distance_matrix, radius=self.bond_lengths[('Li', 'Si')])
        self.cluster_neighbors = self.get_neighbors(distance_matrix = self.cluster_distance_matrix, radius=self.bond_lengths[('Si', 'Si')])

        clusters = self.find_clusters()
        clusters.sort(key=len)
        self.clusters = clusters
        #avg_distance = self.get_mean_distance(clusters)

        return clusters

    def get_distance_matrix(self, structure):
        distance_matrix = [[0 for x in structure.sites] for y in structure.sites]
        for i in range(len(structure.sites)):
            for j in range(i):
                _distance = structure.sites[i].distance(structure.sites[j])
                distance_matrix[j][i] = _distance
                distance_matrix[i][j] = _distance

        return distance_matrix

    def get_track_distance_matrix(self, structure, track_els):
        track_positions=[]
        cluster_positions=[]
        for i in range(len(structure.species)):
            if structure.species[i] in track_els:
                track_positions.append(i)
            else:
                cluster_positions.append(i)

        distance_matrix = [[0 for x in cluster_positions] for y in track_positions]
        for i in range(len(track_positions)):
            for j in range(len(cluster_positions)):
                _distance = structure.sites[track_positions[i]].distance(structure.sites[cluster_positions[j]])
                distance_matrix[i][j] = _distance
        return distance_matrix

    def get_neighbors(self, distance_matrix, radius):
        neighbors = [[] for x in range(len(distance_matrix))]
        for i in range(len(distance_matrix)):
            for j in range(len(distance_matrix[i])):
                if distance_matrix[i][j] <= radius:
                    neighbors[i].append(j)
        return neighbors

    def get_mean_distance(self, clusters):
        avg_dist = [[] for x in clusters]
        for n in range(len(clusters)):
            _cluster = clusters[n]
            bonds = 0
            total_distance = 0
            for i in _cluster:
                for j in _cluster:
                    _distance = self.cluster_distance_matrix[i][j]
                    if _distance > 0:
                        total_distance += _distance
                        bonds += 1
            avg_dist[n] = 0 if bonds == 0 else total_distance / bonds
        return avg_dist


    def find_clusters(self):
        _sites = set(np.arange(len(self.cluster_neighbors)))
        clusters = []
        while _sites:
            site = _sites.pop()
            _cluster = {site}
            queue = [site]
            while queue:
                _site = queue.pop()
                _neighbors = set(self.cluster_neighbors[_site])
                _neighbors.difference_update(_cluster)
                _sites.difference_update(_neighbors)
                _cluster.update(_neighbors)
                queue.extend(_neighbors)
            clusters.append(_cluster)
        return clusters

    def get_clusters_as_structures(self):
        clusters_sites = []
        for cluster in self.clusters:
            cluster_sites = []
            for i in cluster:
                cluster_sites.append(self.input_structure.sites[i])
            clusters_sites.append(cluster_sites)
        cluster_structs = [Structure.from_sites(cluster_sites) for cluster_sites in clusters_sites]
        return cluster_structs

    def get_rdfs(self):
        rdfs = []
        return rdfs

    def plot_clusterSize(self):
        return