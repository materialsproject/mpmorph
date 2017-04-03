from pymatgen import Site, Structure, Element
import numpy as np

class ClusteringAnalyzer(object):
    """

    """

    def __init__(self, radius = 2.6):
        self.radius = 2.6

    def test(self, input_structure, radius):
        '''
        Add founctionality to remove elements other than the desired
        :param input_structure:
        :param radius:
        :return:
        '''
        self.radius = radius
        pruned_structure = input_structure
        pruned_structure.remove_species([Element("Li")])
        struct_sites = pruned_structure.sites
        distance_matrix = self.get_distance_matrix(struct_sites)
        neighbors = self.get_neighbors(distance_matrix, struct_sites, input_structure)
        clusters = self.get_clusters(neighbors)
        clusters.sort(key=len)
        avg_distance = self.get_mean_distance(distance_matrix, clusters)
        return clusters, avg_distance

    def get_distance_matrix(self, sites):
        distance_matrix = [[0 for x in sites] for y in sites]
        for i in range(len(sites)):
            for j in range(i):
                #_distance = np.sqrt((sites[i].x-sites[j].x)**2 + (sites[i].y-sites[j].y)**2 + (sites[i].z-sites[j].z)**2)
                #if _distance >= np.min(structure.lattice.abc)/2:
                _distance = sites[i].distance(sites[j])
                distance_matrix[j][i] = _distance
                distance_matrix[i][j] = _distance

        return distance_matrix

    def get_neighbors(self, distance_matrix, sites, structure):
        neighbors = [[] for x in sites]
        for i in range(len(sites)):
            for j in range(len(sites)):
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


    def get_clusters(self, neighbors):
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