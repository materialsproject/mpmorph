from pymatgen import Site, Structure, Element
import numpy as np

class ClusteringAnalyzer(object):
    """

    """

    def __init__(self, input_structure, radius = 2.6):
        self.radius = 2.6
        self.input_structure = input_structure
        self.clusters = []
        self.struct_sites = []

    def get_clusters_as_structures(self):
        clusters_sites = []
        for cluster in self.clusters:
            cluster_sites = []
            for i in cluster:
                cluster_sites.append(self.struct_sites[i])
            clusters_sites.append(cluster_sites)
        cluster_structs = [ Structure.from_sites(cluster_sites) for cluster_sites in clusters_sites]
        return cluster_structs

    def get_rdfs(self):
        rdfs = []
        return rdfs

    def get_clusters(self, input_structure=None, radius = None):
        '''
        TODO: Add functionality to remove elements other than the desired
              RDF's of Cluster
              Return Sites in Cluster
        :param input_structure:
        :param radius:
        :return:
        '''
        self.radius = radius or self.radius
        self.input_structure = input_structure or self.input_structure

        pruned_structure = input_structure
        pruned_structure.remove_species([Element("Li")])
        self.struct_sites = pruned_structure.sites
        distance_matrix = self.get_distance_matrix()
        neighbors = self.get_neighbors(distance_matrix, input_structure)
        clusters = self.find_clusters(neighbors)
        clusters.sort(key=len)
        self.clusters = clusters
        avg_distance = self.get_mean_distance(distance_matrix, clusters)
        return clusters, avg_distance

    def get_distance_matrix(self):
        distance_matrix = [[0 for x in self.struct_sites] for y in self.struct_sites]
        for i in range(len(self.struct_sites)):
            for j in range(i):
                #_distance = np.sqrt((sites[i].x-sites[j].x)**2 + (sites[i].y-sites[j].y)**2 + (sites[i].z-sites[j].z)**2)
                #if _distance >= np.min(structure.lattice.abc)/2:
                _distance = self.struct_sites[i].distance(self.struct_sites[j])
                distance_matrix[j][i] = _distance
                distance_matrix[i][j] = _distance

        return distance_matrix

    def get_neighbors(self, distance_matrix, structure):
        neighbors = [[] for x in self.struct_sites]
        for i in range(len(self.struct_sites)):
            for j in range(len(self.struct_sites)):
                if distance_matrix[i][j] <= self.radius:
                    neighbors[i].append(j)
                #elif distance_matrix[i][j] >= np.min(structure.lattice.abc)/2:
                    #if sites[i].distance(sites[j]) <= self.radius:
                        #neighbors[i].append((j))
        return neighbors

    def get_mean_distance(self, distance_matrix, clusters):
        avg_dist = [[] for x in clusters]
        for n in range(len(clusters)):
            _cluster = clusters[n]
            bonds = 0
            total_distance = 0
            for i in _cluster:
                for j in _cluster:
                    _distance = distance_matrix[i][j]
                    if _distance > 0:
                        total_distance += _distance
                        bonds += 1
            avg_dist[n] = 0 if bonds == 0 else total_distance / bonds
        return avg_dist


    def find_clusters(self, neighbors):
        _sites = set(np.arange(len(neighbors)))
        clusters = []
        while _sites:
            site = _sites.pop()
            _cluster = {site}
            queue = [site]
            while queue:
                _site = queue.pop()
                _neighbors = set(neighbors[_site])
                _neighbors.difference_update(_cluster)
                _sites.difference_update(_neighbors)
                _cluster.update(_neighbors)
                queue.extend(_neighbors)
            clusters.append(_cluster)
        return clusters

    def plot_clusterSize(self):
        return