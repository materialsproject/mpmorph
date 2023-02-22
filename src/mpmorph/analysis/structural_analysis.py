import itertools
from copy import deepcopy
from multiprocessing import Pool

import numpy as np
from pymatgen.analysis import structure_analyzer
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.util.coord import get_angle
from scipy.spatial import Voronoi


def polyhedra_connectivity(structures, pair, cutoff, step_freq=1):
    """
    Args:
        structures:
        pair:
        cutoff:
        step_freq:
        given: Given polyhedra are of this

    Returns:

    """
    n_frames = len(structures)
    center_atom = pair[0]
    shell_atom = pair[1]

    connectivity = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    connectivity_template = deepcopy(connectivity)
    connectivity_sub_categories = {}

    for s_index in itertools.count(0, step_freq):
        if s_index >= n_frames:
            break
        structure = structures[s_index]

        polyhedra_list = []

        for i in range(len(structure)):
            current_poly = []
            if str(structure[i].specie) == center_atom:
                for j in range(len(structure)):
                    if str(structure[j].specie) == shell_atom:
                        d = structure.get_distance(i, j)
                        if d < cutoff:
                            current_poly.append(j)
            polyhedra_list.append(set(current_poly))

        for polypair in itertools.combinations(polyhedra_list, 2):
            polyhedra_pair_type = (len(polypair[0]), len(polypair[1]))

            shared_vertices = len(polypair[0].intersection(polypair[1]))

            if shared_vertices in connectivity:
                connectivity[shared_vertices] += 1

            if shared_vertices:
                if polyhedra_pair_type in connectivity_sub_categories:
                    if (
                        shared_vertices
                        in connectivity_sub_categories[polyhedra_pair_type]
                    ):
                        connectivity_sub_categories[polyhedra_pair_type][
                            shared_vertices
                        ] += 1
                elif polyhedra_pair_type[::-1] in connectivity_sub_categories:
                    if (
                        shared_vertices
                        in connectivity_sub_categories[polyhedra_pair_type[::-1]]
                    ):
                        connectivity_sub_categories[polyhedra_pair_type[::-1]][
                            shared_vertices
                        ] += 1
                else:
                    connectivity_sub_categories[polyhedra_pair_type] = deepcopy(
                        connectivity_template
                    )
                    if (
                        shared_vertices
                        in connectivity_sub_categories[polyhedra_pair_type]
                    ):
                        connectivity_sub_categories[polyhedra_pair_type][
                            shared_vertices
                        ] = 1
    return connectivity, connectivity_sub_categories


def coordination_number_distribution(structures, pair, cutoff, step_freq=1):
    """
    Calculates coordination number distribution
    Args:
        structures:
        pair:
        cutoff:
        step_freq:

    Returns:

    """

    cn_list = []
    n_frames = len(structures)
    for s_index in itertools.count(0, step_freq):
        if s_index >= n_frames:
            break
        structure = structures[s_index]
        for i in range(len(structure)):
            if str(structure[i].specie) == pair[0]:
                cn = 0
                for j in range(len(structure)):
                    if str(structure[j].specie) == pair[1]:
                        d = structure.get_distance(i, j)
                        if d < cutoff and np.abs(d) > 0.1:
                            cn += 1
                cn_list.append(cn)
    return cn_list


def get_cn(structure, pair, cutoff):
    cn_list = []
    for i in range(len(structure)):
        if str(structure[i].specie) == pair[0]:
            cn = 0
            for j in range(len(structure)):
                if str(structure[j].specie) == pair[1]:
                    d = structure.get_distance(i, j)
                    if d < cutoff and np.abs(d) > 0.1:
                        cn += 1
            cn_list.append(cn)
    return cn_list


class BondAngleDistribution(object):
    """
    Bond Angle Distribution
    Args:
        structures (list): list of structures
        cutoffs (dict): a dictionary of cutoffs where keys are tuples of pairs ('A','B')
        step_freq: calculate every this many steps
    Attributes:
        bond_angle_distribution (dict)
    """

    def __init__(self, structures, cutoffs, step_freq=1):
        self.bond_angle_distribution = None
        self.structures = structures
        self.step_freq = step_freq
        self.unique_triplets = self.get_unique_triplets(structures[0])
        if isinstance(cutoffs, dict):
            self.cutoffs = cutoffs
            self._cutoff_type = "dict"
        elif isinstance(cutoffs, float):
            self.cutoffs = cutoffs
            self._cutoff_type = "constant"
        else:
            raise ValueError(
                "Cutoffs must be specified as dict of pairs or globally as a single flaot."
            )

    @property
    def n_frames(self):
        return len(self.structures)

    def get_angle(self, s_index, i, j, k):
        """
        Returns **Minimum Image** angle specified by three sites.

        Args:
            s_index: Structure index in structures list
            i (int): Index of first site.
            j (int): Index of second site.
            k (int): Index of third site.

        Returns:
            (float) Angle in degrees.
        """
        structure = self.structures[s_index]
        lat_vec = np.array(
            [structure.lattice.a, structure.lattice.b, structure.lattice.c]
        )

        v1 = structure[i].coords - structure[j].coords
        v2 = structure[k].coords - structure[j].coords

        for v in range(3):
            if np.fabs(v1[v]) > lat_vec[v] / 2.0:
                v1[v] -= np.sign(v1[v]) * lat_vec[v]
            if np.fabs(v2[v]) > lat_vec[v] / 2.0:
                v2[v] -= np.sign(v2[v]) * lat_vec[v]
        return get_angle(v1, v2, units="degrees")

    @staticmethod
    def get_unique_triplets(s):
        central_atoms = s.symbol_set
        import itertools

        possible_end_members = []
        for i in itertools.combinations_with_replacement(central_atoms, 2):
            possible_end_members.append(i)
        unique_triplets = []
        for i in central_atoms:
            for j in possible_end_members:
                triplet = (j[0], i, j[1])
                unique_triplets.append(triplet)
        return unique_triplets

    def _check_skip_triplet(self, s_index, i, n1, n2):
        """
        Helper method to find if a triplet should be skipped
        Args:
            s_index: index of structure in self.structures
            i: index of the central site
            n1: index of the first neighbor site
            n2: index of the second neighbor site
        Returns:
            True if pair distance is longer than specified in self.cutoffs
        """
        ns = [n1, n2]
        s = self.structures[s_index]
        skip_triplet = False
        for j in ns:
            pair = (s[i].species_string, s[j].species_string)
            if pair not in self.cutoffs:
                pair = pair[::-1]
            if s.get_distance(i, j) > self.cutoffs[pair]:
                skip_triplet = True
                break
        return skip_triplet

    def get_bond_angle_distribution(self):
        bond_angle_dict = {}
        for triplet in self.unique_triplets:
            bond_angle_dict[triplet] = np.zeros(180 + 1)

        for s_index in itertools.count(0, self.step_freq):
            if s_index >= self.n_frames:
                break
            s = self.structures[s_index]

            # Narrow down the search space around a given atom for neighbors
            if self._cutoff_type == "dict":
                neighbor_search_cutoff = max(self.cutoffs.values())
            else:
                neighbor_search_cutoff = self.cutoffs
            neighbors = s.get_all_neighbors(neighbor_search_cutoff, include_index=True)

            for i in range(len(s)):
                el_origin = s[i].species_string

                # get all pair combinations of neoghbor sites of i:
                for p in itertools.combinations(neighbors[i], 2):
                    # check if pairs are within the defined cutoffs
                    if self._cutoff_type == "dict":
                        if self._check_skip_triplet(s_index, i, p[0][2], p[1][2]):
                            continue
                    else:
                        pass

                    angle = self.get_angle(s_index, p[0][2], i, p[1][2])

                    # round to nearest integer
                    angle = int(np.rint([angle])[0])
                    el1 = p[0][0].species_string
                    el2 = p[1][0].species_string
                    if (el1, el_origin, el2) in bond_angle_dict:
                        bond_angle_dict[(el1, el_origin, el2)][angle] += 1
                    elif (el2, el_origin, el1) in bond_angle_dict:
                        bond_angle_dict[(el2, el_origin, el1)][angle] += 1
                    else:
                        print(el1, el_origin, el2)
                        raise KeyError("Problem finding the triplet!!!")
        for triplet in bond_angle_dict:
            total = np.sum(bond_angle_dict[triplet])
            if total != 0.0:
                bond_angle_dict[triplet] /= total
        self.bond_angle_distribution = bond_angle_dict

    def plot_bond_angle_distribution(self):
        import matplotlib.pyplot as plt

        if not self.bond_angle_distribution:
            self.get_bond_angle_distribution()
        plt.figure()
        triplets = self.bond_angle_distribution.keys()
        legend = []
        for trip in triplets:
            legend.append("-".join(trip))
        for triplet in triplets:
            plt.plot(range(180 + 1), self.bond_angle_distribution[triplet])

        plt.xlabel("Angle (degrees)")
        plt.ylabel("Frequency (fractional)")
        plt.legend(legend, loc=0)
        return plt

    def get_binary_angle_dist_plot(self, title=None):
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 6))
        c = 0
        maxes = []
        for triplet in self.bond_angle_distribution:
            p = self.bond_angle_distribution[triplet]
            c += 1
            ax = fig.add_subplot(2, 3, c)
            ax.plot(range(len(p)), p)
            maxes.append(max(p))
            ax.annotate(
                "-".join(triplet), (0.75, 0.88), xycoords="axes fraction", size=16
            )
            ax.set_yticklabels([])
            ax.xaxis.set_ticks(np.arange(0, 181, 30))
            plt.gca().set_ylim([0, 0.1])
            if c in [1, 2, 3]:
                ax.set_xticklabels([])
            else:
                plt.xlabel("Angle (degrees)", fontsize=16)
            if c in [1, 4]:
                plt.ylabel("Intensity (a.u.)", fontsize=16)
            ax.tick_params(axis="both", which="major", labelsize=16)
            if c == 2:
                if title:
                    plt.title(title)
        for ax in fig.axes:
            ax.set_ylim(0.0, max(maxes))
        fig.subplots_adjust(top=0.75)
        fig.tight_layout()
        return fig


def compute_mean_coord(structures, freq=100):
    """
    NOTE: This function will be removed as it has been migrated
    to pymatgen.
    Calculate average coordination numbers
    With voronoi polyhedra construction
    args:
        - structures: list of Structures
        - freq: sampling frequency of coord number [every freq steps]
    returns:
        - a dictionary of elements and corresponding mean coord numbers
    """
    cn_dict = {}
    for el in structures[0].composition.elements:
        cn_dict[el.name] = 0.0
    count = 0
    for t in range(len(structures)):
        if t % freq != 0:
            continue
        count += 1
        vor = structure_analyzer.VoronoiCoordFinder(structures[t])
        for atom in range(len(structures[0])):
            CN = vor.get_coordination_number(atom)
            cn_dict[structures[t][atom].species_string] += CN
    el_dict = structures[0].composition.as_dict()
    for el in cn_dict:
        cn_dict[el] = cn_dict[el] / el_dict[el] / count
    return cn_dict


class VoronoiAnalysis(object):
    """
    NOTE: This class has also been migrated to pymatgen so will be removed!
    """

    def __init__(self):
        self.vor_ensemble = None

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
        vor_index = np.array([0, 0, 0, 0, 0, 0, 0, 0])

        for key in voro.ridge_dict:
            if 0 in key:
                "This means if the center atom is in key"
                if -1 in key:
                    "This means if an infinity point is in key"
                    print("Cutoff too short. Exiting.")
                    return None
                else:
                    try:
                        vor_index[len(voro.ridge_dict[key]) - 3] += 1
                    except IndexError:
                        # If a facet has more than 10 edges, it's skipped here
                        pass
        return vor_index

    def from_structures(
        self, structures, cutoff=4.0, step_freq=10, qhull_options="Qbb Qc Qz"
    ):
        """
        A constructor to perform Voronoi analysis on a list of pymatgen structrue objects

        Args:
            structures (list): list of Structures
            cutoff (float: cutoff distance around an atom to search for neighbors
            step_freq (int): perform analysis every step_freq steps
            qhull_options (str): options to pass to qhull
        Returns:
            A list of [voronoi_tesellation, count]
        """
        print("This might take a while...")
        voro_dict = {}
        step = 0
        for structure in structures:
            step += 1
            if step % step_freq != 0:
                continue

            v = []
            for n in range(len(structure)):
                v.append(
                    str(
                        self.voronoi_analysis(
                            structure, n=n, cutoff=cutoff, qhull_options=qhull_options
                        ).view()
                    )
                )
            for voro in v:
                if voro in voro_dict:
                    voro_dict[voro] += 1
                else:
                    voro_dict[voro] = 1
        self.vor_ensemble = sorted(
            voro_dict.items(), key=lambda x: (x[1], x[0]), reverse=True
        )[:15]
        return self.vor_ensemble

    @property
    def plot_vor_analysis(self):
        import matplotlib.pyplot as plt

        t = zip(*self.vor_ensemble)
        labels = t[0]
        val = list(t[1])
        tot = np.sum(val)
        val = [float(j) / tot for j in val]
        pos = np.arange(len(val)) + 0.5  # the bar centers on the y axis
        plt.figure(figsize=(4, 4))
        plt.barh(pos, val, align="center", alpha=0.5)
        plt.yticks(pos, labels)
        plt.xlabel("Fraction")
        plt.title("Voronoi Spectra")
        plt.grid(True)
        return plt


class RadialDistributionFunction(object):
    """
    Class to calculate partial radial distribution functions (RDFs) of sites.
    Typically used to analyze the pair correlations in liquids or amorphous structures.
    Supports multiprocessing: see get_radial_distribution_function method

    Args:
        structures (list): a list of Structure objects.
        cutoff (float): maximum distance to search for pairs. (defauly = 5.0)
            Note cutoff should be smaller than or equal to the half of the edge
            length of the box due to periodic boundaires.
        bin_size (float): thickness of each coordination shell in Angstroms (default = 0.1)
        step_freq (int): compute and store RDFs every step_freq steps
            to average later. (default = 2)
        smooth (int): number of smoothing passes (default = 1)
        title (str): title for the RDF plot.
    Returns:
        A dictionary of partial radial distribution functions with pairs as keys and RDFs as values.
        RDFs themselves are arrays of length cutoff/bin_size.
    """

    def __init__(
        self,
        structures,
        cutoff=5.0,
        bin_size=0.1,
        step_freq=2,
        smooth=1,
        title="Radial distribution functions",
    ):
        self.structures = structures
        self.cutoff = cutoff
        self.bin_size = bin_size
        self.step_freq = int(step_freq)
        self.smooth = smooth
        self.n_frames = int(len(self.structures))
        self.n_atoms = len(self.structures[0])
        self.n_species = self.structures[0].composition.as_dict()
        self.get_pair_order = None
        self.title = title
        self.RDFs = {}
        ss = self.structures[0].symbol_set
        self.pairs = [p for p in itertools.combinations_with_replacement(ss, 2)]
        self.counter = 1

    @property
    def n_bins(self):
        _bins = int(np.ceil(self.cutoff / self.bin_size))
        if _bins < 2:
            raise ValueError("More bins required!")
        return _bins

    def get_radial_distribution_functions(self, nproc=1):
        """
        Args:
            nproc: (int) number of processors to utilize (defaults to 1)
        Returns:
            A dictionary of partial radial distribution functions
            with pairs as keys and RDFs as values.
            Each RDF arrays of length cutoff/bin_size.
        """

        frames = [
            (
                self.structures[i * self.step_freq],
                self.pairs,
                self.n_bins,
                self.cutoff,
                self.bin_size,
            )
            for i in range(int(self.n_frames / self.step_freq))
        ]
        self.counter = len(frames)
        pool = Pool(nproc)
        results = pool.map(_process_frame, frames)
        pool.close()
        pool.join()

        # Collect all rdfs
        for pair in self.pairs:
            self.RDFs[pair] = np.zeros(self.n_bins)
        for i in results:
            for pair in self.pairs:
                self.RDFs[pair] += i[pair]

        self.get_pair_order = []

        for i in self.RDFs.keys():
            self.get_pair_order.append("-".join(list(i)))
            density_of_atom2 = self.n_species[i[1]] / self.structures[0].volume
            for j in range(self.n_bins):
                r = j * self.bin_size
                if r == 0:
                    continue
                self.RDFs[i][j] = (
                    self.RDFs[i][j]
                    / self.n_species[i[0]]
                    / 4.0
                    / np.pi
                    / r
                    / r
                    / self.bin_size
                    / density_of_atom2
                    / self.counter
                )

        if self.smooth:
            self.RDFs = get_smooth_rdfs(self.RDFs, passes=self.smooth)
        return self.RDFs

    def plot_radial_distribution_functions(self):
        """
        :return: a plot of RDFs
        """
        import matplotlib.pyplot as plt

        x = []
        for j in range(self.n_bins):
            r = j * self.bin_size
            x.append(r)
        rdfs = self.RDFs
        plt.figure()
        for rdf in rdfs:
            plt.plot(x, rdfs[rdf])

        plt.xlabel("$r$, distance (Angstrom)")
        plt.ylabel("g($r$)")
        plt.legend(
            self.get_pair_order,
            bbox_to_anchor=(0.975, 0.975),
            loc=0,
            borderaxespad=0.0,
            prop={"family": "sans-serif", "size": 13},
        )
        plt.title(self.title)
        return plt


def _process_frame(data):
    """
    Helper function for parallel rdf computation
    """
    coord_frame, pairs, n_bins, cutoff, bin_size = (
        data[0],
        data[1],
        data[2],
        data[3],
        data[4],
    )
    process_RDFs = {}
    n_atoms = len(coord_frame)

    for pair in pairs:
        process_RDFs[pair] = np.zeros(n_bins)

    distance_matrix = coord_frame.distance_matrix

    for atom1 in range(n_atoms - 1):
        atom1_specie = coord_frame[atom1].species_string
        for atom2 in range(atom1 + 1, n_atoms):
            atom2_specie = coord_frame[atom2].species_string
            if distance_matrix[atom1, atom2] > cutoff:
                continue
            bin_index = int(distance_matrix[atom1, atom2] / bin_size)
            # skip itself
            if bin_index == 0:
                continue
            key = (atom1_specie, atom2_specie)
            if key in process_RDFs:
                process_RDFs[key][bin_index] += 1
            if key[::-1] in process_RDFs:
                process_RDFs[key[::-1]][bin_index] += 1
    return process_RDFs


def get_smooth_rdfs(RDFs, passes=1):
    """
    Helper function to recursively smooth RDFs using a 5-parameter Savitzky-Golay filter.
    Args:
        RDFs: A dictionary of partial radial distribution functions
        with pairs as keys and RDFs as values.
        passes: number of times the filter is applied during smoothing.
    Returns
        RDFs dictionary with with each RDF smoothed.
    """
    if passes == 0:
        return RDFs
    else:
        for rdf in RDFs:
            smooth_RDF = deepcopy(RDFs[rdf])
            for j in range(2, len(RDFs[rdf]) - 2):
                smooth_RDF[j] = (
                    -3 * RDFs[rdf][j - 2]
                    + 12 * RDFs[rdf][j - 1]
                    + 17 * RDFs[rdf][j]
                    + 12 * RDFs[rdf][j + 1]
                    - 3 * RDFs[rdf][j + 2]
                ) / 35.0
            RDFs[rdf] = smooth_RDF
        passes -= 1
        return get_smooth_rdfs(RDFs, passes=passes)


def get_smooth_rdfs(RDFs, passes=1):
    """
    Helper function to recursively smooth RDFs using a 5-parameter Savitzky-Golay filter.
    Args:
        RDFs: A dictionary of partial radial distribution functions
        with pairs as keys and RDFs as values.
        passes: number of times the filter is applied during smoothing.
    Returns
        RDFs dictionary with with each RDF smoothed.
    """
    if passes == 0:
        return RDFs
    else:
        for rdf in RDFs:
            smooth_RDF = deepcopy(RDFs[rdf])
            for j in range(2, len(RDFs[rdf]) - 2):
                smooth_RDF[j] = (
                    -3 * RDFs[rdf][j - 2]
                    + 12 * RDFs[rdf][j - 1]
                    + 17 * RDFs[rdf][j]
                    + 12 * RDFs[rdf][j + 1]
                    - 3 * RDFs[rdf][j + 2]
                ) / 35.0
            RDFs[rdf] = smooth_RDF
        passes -= 1
        return get_smooth_rdfs(RDFs, passes=passes)


def get_sample_structures(structures, n=10, steps_skip_first=1000):
    """
    Helper method to extract n unique structures from an MD output
    Args:
        xdatcar_path: (str) path
    Returns:
        A list of Structure objects
    """
    if type(structures) == type("asdf"):
        input_structures = Xdatcar(xdatcar_path).structures
    else:
        input_structures = structures
    output_structures = []
    t = len(input_structures) - steps_skip_first
    for i in range(n):
        output_structures.append(input_structures[::-1][i * t // n])
    return output_structures
