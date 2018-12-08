'''
Created on Feb 13, 2014

@author: sushant
'''


class Cluster(object):
    def __init__(self,  name, dim):
        self.name = name
        self.dim = dim
        self._points = []

    @property
    def points(self):
        return self._points

    def add_point(self, point):
        self._points.append(point)

    def get_x(self):
        return [p[0] for p in self._points]
    
    def get_y(self):
        return [p[1] for p in self._points]

    def get_z(self):
        if self.dim > 2:
            return [p[2] for p in self._points]
        return None
    
    def has(self, point):
        return point in self._points
            
    def __str__(self):
        return "%s: %d points" % (self.name, len(self._points))