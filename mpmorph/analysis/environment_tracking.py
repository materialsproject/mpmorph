from pymatgen import Structure, Composition
from pymatgen.io.vasp import Xdatcar
from mpmorph.analysis.clustering_analysis import ClusteringAnalyzer
from mpmorph.analysis.structural_analysis import RadialDistributionFunction

class EnvironmentTracker():

    def __init__(self):
        pass

    def run(self, xdatcar_path):
        xdat = Xdatcar(xdatcar_path)
        structures = xdat.structures
        neighbor_array = [] #Need to predeclare 3d array of length of xdatcar
        for i in range(len(structures)):
            frame_num, neighbors = process_frame(i, structures[i])
        return

    def process_frame(self, frame_num, structure):

        return frame_num, neighbors