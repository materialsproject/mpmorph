from pymatgen.analysis import structure_analyzer
from scipy.spatial import Voronoi
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import itertools

def compute_mean_coord(structures, freq = 100):
    '''
    Calculate average coordination numbers
    With voronoi polyhedra construction
    args:
        - structures: list of pymatgen structures
        - freq: sampling frequency of coord number [every freq steps]
    returns:
        - a dictionary of elements and corresponding mean coord numbers
    '''
    cn_dict={}
    for el in structures[0].composition.elements:
        cn_dict[el.name]=0.0
    count = 0
    for t in range(len(structures)):
        if t%freq != 0:
            continue
        count += 1
        vor = structure_analyzer.VoronoiCoordFinder(structures[t])
        for atom in range(len(structures[0])):
            CN = vor.get_coordination_number(atom)
            cn_dict[structures[t][atom].species_string] += CN
    el_dict = structures[0].composition.as_dict()
    for el in cn_dict:
        cn_dict[el] = cn_dict[el]/el_dict[el]/count
    return cn_dict

class VoronoiAnalysis(object):

    def __init__(self):
        self.vor_ensemble = None
        pass

    @staticmethod
    def voronoi_analysis(structure, n=0, cutoff=5.0, qhull_options="Qbb Qc Qz"):
        """
        Performs voronoi analysis and returns the polyhedra around atom n
        in Schlaefli notation.

        Note: Part of this function is obtained from the pymatgen VoronoiCoordinationFinder
              and should be merged with that function in future.

        Args:
            - structure: pymatgen Structure object
            - n: index of the center atom in structure
            - cutoff: cutoff distance around n to search for neighbors
        Returns:
            - voronoi index of n: <x3,x4,x5,x6,x6,x7,x8,x9,x10>
              where x_i denotes number of facets with i edges
        """

        center = structure[n]
        neighbors = structure.get_sites_in_sphere(center.coords, cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = np.array([s.coords for s in neighbors])
        voro = Voronoi(qvoronoi_input, qhull_options=qhull_options)
        vor_index = np.array([0,0,0,0,0,0,0,0])

        for key in voro.ridge_dict:
            if 0 in key:
                "This means if the center atom is in key"
                if -1 in key:
                    "This means if an infinity point is in key"
                    print "Cutoff too short. Exiting."
                    return None
                else:
                    try:
                        vor_index[len(voro.ridge_dict[key])-3]+=1
                    except IndexError:
                        # If a facet has more than 10 edges, it's skipped here
                        pass
        return vor_index



    def from_structures(self, structures, cutoff=4.0, step_freq=10):
        """
        A constructor to perform Voronoi analysis on a list of pymatgen structrue objects
        """
        print "This might take a while..."
        voro_dict = {}
        step = 0
        for structure in structures:
            step+=1
            if step%step_freq != 0:
                continue

            v = []
            for n in range(len(structure)):
                v.append(str(self.voronoi_analysis(structure,n=n,cutoff=cutoff, qhull_options="Qbb Qc Qz").view()))
            for voro in v:
                if voro in voro_dict:
                    voro_dict[voro]+=1
                else:
                    voro_dict[voro]=1
        self.vor_ensemble = sorted(voro_dict.items(), key=lambda x: (x[1],x[0]), reverse=True)[:15 ]
        return self.vor_ensemble

    @property
    def plot_vor_analysis(self):
        t = zip(*self.vor_ensemble)
        labels = t[0]
        val = list(t[1])
        tot = np.sum(val)
        val = [float(j)/tot for j in val]
        pos = np.arange(len(val))+.5    # the bar centers on the y axis
        plt.figure()
        plt.barh(pos,val, align='center', alpha=0.5)
        plt.yticks(pos, labels)
        plt.xlabel('Count')
        plt.title('Voronoi Spectra')
        plt.grid(True)
        return plt

class RDF(object):
    # RDF
    # bin_size = 0.1 Angstrom, default
    # cutoff = 10 Angstrom, default
    # step_freq = 2, default
    # smooth = # of passes for a 5-point Savitzky-Golay signal smoothing
    # Find unique species and generate pairs
    # For Cr-O system, Cr-Cr, O-O, Cr-O, and O-Cr
    # Last two RDFs are similar, and they are both
    # calculated to allow the user to select the desired one


    def __init__(self, xdatcar, cutoff = 5.0, bin_size = 0.025, step_freq = 2, smooth = 1, title="Radial distribution functions\n(T = 1400 K)"):
        '''
        :param xdatcar:
        :param cutoff:
        :param bin_size:
        :param step_freq:
        :param smooth:
        :return:
        '''
        self.xdatcar = xdatcar
        self.cutoff = cutoff
        self.bin_size = bin_size
        self.step_freq = step_freq
        self.smooth = smooth
        self.n_frames = len(xdatcar.structures)
        self.n_atoms = len(xdatcar.structures[0])
        self.n_species = xdatcar.structures[0].composition.as_dict()
        self.get_pair_order = None
        self.title = title

    @property
    def n_bins(self):
        _bins = int(self.cutoff/self.bin_size)
        if _bins <2:
            raise ValueError("More bins required!")
        return _bins

    def compute_RDF(self):
        '''
        :return:
        '''
        self.RDFs = {}

        ss = self.xdatcar.structures[0].symbol_set
        self.pairs = itertools.combinations_with_replacement(ss,2)
        for pair in self.pairs:
            self.RDFs[pair]=np.zeros(self.n_bins)

        counter = 0
        for frame in itertools.count(0, self.step_freq):
            if frame >= self.n_frames:
                break
            counter +=1
            # Coordinates in the current frame
            coord_frame = self.xdatcar.structures[frame]
            distance_matrix = coord_frame.distance_matrix
            for atom1 in range(self.n_atoms):
                for atom2 in range(self.n_atoms):
                    atom1_specie = coord_frame[atom1].species_string
                    atom2_specie = coord_frame[atom2].species_string
                    if distance_matrix[atom1,atom2] > self.cutoff:
                        continue
                    bin_index = int(distance_matrix[atom1,atom2]/self.bin_size)
                    #skip itself
                    if bin_index == 0:
                        continue
                    key = (atom1_specie,atom2_specie)
                    if key in self.RDFs:
                        self.RDFs[key][bin_index] += 1
        self.totalRDFs = deepcopy(self.RDFs)
        self.get_pair_order = []
        for i in self.RDFs.keys():
            self.get_pair_order.append('-'.join(list(i)))
            density_of_atom2 = self.n_species[i[1]]/self.xdatcar.structures[0].volume
            for j in range(self.n_bins):
                r = j*self.bin_size
                if r == 0:
                    continue
                # Divide by number of atom1 to obtain the average of atom2 at r+dr per atom type 1
                # Divide by 4pi r^2 dr to convert to density
                # Divide by counter to get average in time
                self.totalRDFs[i][j] = self.totalRDFs[i][j]/self.n_species[i[0]]/counter
                self.RDFs[i][j]=self.RDFs[i][j]/self.n_species[i[0]]/4/np.pi/r/r/self.bin_size/density_of_atom2/counter

        if self.smooth:
            self.RDFs= get_smooth_RDF(self.RDFs,passes=self.smooth)
        return self.RDFs

    def plot_RDF(self,total=False):
        x = []
        for j in range(self.n_bins):
            r = j*self.bin_size
            x.append(r)
        if total:
            rdfs = self.totalRDFs
        else:
            rdfs = self.RDFs
        for rdf in rdfs:
            plt.plot(x,rdfs[rdf])

        plt.xlabel("$r$, distance (Angstrom)")
        plt.ylabel("g($r$)")
        plt.legend(self.get_pair_order,bbox_to_anchor=(0.975, 0.975), loc=1,
                   borderaxespad=0.,prop={'family': 'sans-serif', 'size':13})
        plt.title(self.title)
        return plt


def get_smooth_RDF(RDFs,passes=1):
    #Smooth?
    """
    :param passes:
    :return:
    """
    if passes==0:
       return RDFs
    else:
        print 'Smoothing. pass remaining: ', passes
        for rdf in RDFs:
            smooth_RDF = deepcopy(RDFs[rdf])
            for j in range(2,len(RDFs[rdf])-2):
                smooth_RDF[j] = (-3*RDFs[rdf][j-2]+12*RDFs[rdf][j-1]
                                 +17*RDFs[rdf][j]+12*RDFs[rdf][j-1]-3*RDFs[rdf][j-2])/35.0
            RDFs[rdf] = smooth_RDF
        passes-=1
        return get_smooth_RDF(RDFs,passes=passes)