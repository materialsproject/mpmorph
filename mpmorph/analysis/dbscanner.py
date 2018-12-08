'''
Created on Feb 13, 2014

@author: sushant
'''

import csv
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpmorph.analysis.cluster import Cluster


class DBScanner:

    def __init__(self, data, eps=0.2, min_pts=2, dim=3):
        self.data = data
        self.eps = eps
        self.min_pts = min_pts
        self.dim = dim
        self.clusters = []
        self.cluster_count = 0
        self.visited = []

        # default noise cluster
        self.noise = Cluster('noise', self.dim)
        self.clusters.append(self.noise)

        for point in self.data:
            if point not in self.visited:
                self.visited.append(point)
                neighbour_pts = self.region_query(point)
                if len(neighbour_pts) < self.min_pts:
                    self.noise.add_point(point)
                else:
                    name = 'cluster-%d' % self.cluster_count
                    new_cluster = Cluster(name, self.dim)

                    self.cluster_count += 1
                    self.expand_cluster(new_cluster, point, neighbour_pts)

    @staticmethod
    def from_file(filename, eps=0.2, min_pts=2, dim=3):
        data = []
        with open(filename, 'r') as file_obj:
            csv_reader = csv.reader(file_obj)
            for id_, row in enumerate(csv_reader):
                point = {'id': id_}
                for i in range(len(row)):
                    point[i] = float(row[i])
                data.append(point)
        return DBScanner(data, eps, min_pts, dim)

    def get_plot(self):
        # Setting up the plot
        fig = plt.figure()

        axis_proj = 'rectilinear'
        if self.dim > 2:
            axis_proj = '%dd' % self.dim

        ax = fig.add_subplot(111, projection=axis_proj)

        for cluster in self.clusters:
            if cluster.name == 'noise':
                continue
            if self.dim == 2:
                ax.scatter(cluster.get_x(), cluster.get_y(), marker='o')
            elif self.dim == 3:
                ax.scatter(cluster.get_x(), cluster.get_y(),
                           cluster.get_z(), marker='o')
        if len(self.noise.points) != 0:
            if self.dim > 2:
                ax.scatter(self.noise.get_x(), self.noise.get_y(), self.noise.get_z(),
                           marker='x', label=self.noise.name)
            else:
                ax.scatter(self.noise.get_x(), self.noise.get_y(),
                           marker='x', label=self.noise.name)

        print("Number of clusters found: %d" % self.cluster_count)
        ax.grid(True)
        plt.show()

    def expand_cluster(self, cluster, point, neighbour_pts):
        cluster.add_point(point)
        for p in neighbour_pts:
            if p not in self.visited:
                self.visited.append(p)
                np = self.region_query(p)
                if len(np) >= self.min_pts:
                    for n in np:
                        if n not in neighbour_pts:
                            neighbour_pts.append(n)

                for other_cluster in self.clusters:
                    if not other_cluster.has(p):
                        if not cluster.has(p):
                            cluster.add_point(p)

                if self.cluster_count == 0:
                    if not cluster.has(p):
                        cluster.add_point(p)

        self.clusters.append(cluster)

    def get_distance(self, from_point, to_point):
        p1 = [from_point[k] for k in range(self.dim)]
        p2 = [to_point[k] for k in range(self.dim)]
        return distance.euclidean(p1, p2)

    def region_query(self, point):
        result = []
        for d_point in self.data:
            if d_point != point:
                if self.get_distance(d_point, point) <= self.eps:
                    result.append(d_point)
        return result
