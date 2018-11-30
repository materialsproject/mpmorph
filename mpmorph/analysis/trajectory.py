from copy import deepcopy
import numpy as np
import bisect
import os
from monty.json import MSONable
from monty.functools import lru_cache
from monty.re import regrep
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen import Structure


class Trajectory(MSONable):
    def __init__(self, base_structure, disp, lattices, frac_coords):
        self.base_structure = base_structure
        self.base_frac_coords = base_structure.frac_coords
        self.disp = np.array(disp)
        self.lattices = lattices
        self.index = 0
        self.structure = base_structure
        self.frac_coords = frac_coords

    def __getattr__(self, attr):

        if hasattr(self.structure, attr):
            return self.structure.__getattr__(attr)
        else:
            raise Exception("Neither Trajectory nor structure has attribute: {}".format(attr))

    def change_index(self, index):
        if index > len(self.disp):
            raise Exception
        else:
            coords = self.base_structure.frac_coords + self.disp[index]
            self.structure = Structure(self.base_structure.lattice,
                                       self.base_structure.species, coords)
            self.index = index

    def change_next(self):
        if self.index + 1 < len(self.disp):
            self.change_index(self.index+1)
        else:
            raise Exception

    def change_previous(self):
        if self.index > 0:
            self.change_index(self.index+1)
        else:
            raise Exception

    def combine(self, trajectory):
        if trajectory.base_structure.lattice != self.base_structure.lattice:
            raise Exception("Lattices are incompatible")
        if trajectory.base_structure.species != self.base_structure.species:
            raise Exception("Elements are not consistent between both trajectories")

        traj_disp = trajectory.disp + self.disp[:, -1:, :]
        self.disp = np.concatenate((self.disp, traj_disp), axis=1)
        if not np.array_equal(self.lattices, trajectory.lattices):
            self.lattices = np.concatenate(self.lattices, trajectory.lattices)

    @lru_cache()
    def as_structures(self):
        structures = [0]*len(self.disp)
        for i in range(len(self.disp)):
            self.change_index(i)
            structures[i] = self.structure.copy()
        return structures

    @property
    @lru_cache()
    def sequential_displacements(self, skip=1):
        seq_displacements = np.subtract(self.disp[::skip],
                                        np.roll(self.disp[::skip], 1, axis=0))
        return seq_displacements

    @classmethod
    def from_structures(cls, structures):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        """
        p, l = [], []
        structure = structures[0]
        for i, s in enumerate(structures):
            p.append(np.array(s.frac_coords)[:, None])
            l.append(s.lattice.matrix)
        p.insert(0, p[0])
        l.insert(0, l[0])
        p = np.concatenate(p, axis=1)  # [site, time step, axis]
        dp = p[:, 1:] - p[:, :-1]
        dp = dp - np.round(dp)
        f_disp = np.cumsum(dp, axis=1)
        c_disp = []
        for i in f_disp:
            c_disp.append([np.dot(d, m) for d, m in zip(i, l[1:])])
        disp = np.array(c_disp)

        # If is NVT-AIMD, clear lattice data.
        if np.array_equal(l[0], l[-1]):
            l = np.array([l[0]])
        else:
            l = np.array(l)

        return cls(structure, disp, l, p[:, 1:])

    @classmethod
    def from_ionic_steps(cls, ionic_steps_dict):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        :param ionic_steps_dict:
        :return:
        """
        return cls.from_structures([Structure.from_dict(i['structure'])
                                    for i in ionic_steps_dict])

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": self.structure.as_dict(),
             "displacements": self.disp.tolist(),
             "lattices": self.lattices.tolist(),
             "frac_coords": self.frac_coords.tolist()
             }
        return d

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, d["displacements"], d["lattices"], d["frac_coords"])


class TemperingTrajectory(MSONable):
    def __init__(self, data, trajectories):
        self.data = data
        self.trajectories = trajectories

        self.index = 0

    def get_structures(self):
        structures = {}
        for key, trajectory in self.trajectories.items():
            structures[key] = trajectory.structure
        return structures

    def change_index(self, index):
        try:
            for key, trajectory in self.trajectories.items():
                trajectory.change_index(index)
            self.index = index
        except:
            raise Exception("Requested index is out of bounds")

    def change_previous(self, index):
        try:
            for key, trajectory in self.trajectories.items():
                trajectory.change_previous(index)
            self.index = index
        except:
            raise Exception("Error: Already at the start of the trajectory")

    def swap_probability(self):
        return self.data["acceptance"][-1]

    def get_trajectory(self, temperature):
        _index = self.get_temperatures().index(temperature)
        return self.trajectories[_index]

    def get_temperatures(self):
        return self.data["temperatures"]

    def plot_swaps(self, show_attempts=False):
        from matplotlib import pyplot as plt
        image_trajectories = self.data["image_trajectories"]
        for i in np.transpose(image_trajectories):
            plt.plot(i)

        if show_attempts:
            for i in self.data["nswap"]:
                plt.plot([i, i], [min(self.get_temperatures()), max(self.get_temperatures())], '--k', linewidth=1)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("time step (#)", fontsize=18)
        plt.ylabel("Temperature", fontsize=18)
        plt.show()

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["data"] = self.data
        d["trajectories"] = [trajectory.as_dict() for trajectory in self.trajectories]
        return d

    @classmethod
    def from_dict(cls, d):
        trajectories = [Trajectory.from_dict(trajectory_dict) for trajectory_dict in d["trajectories"]]
        return cls(d["data"], trajectories)

    @classmethod
    def from_vasp_run(cls, directory, load_structs=False):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        :param ionic_steps_dict:
        :return:
        """

        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        nimages = sum(os.path.isdir(i) and is_number(i) for i in os.listdir(directory))

        raw_data = parse_outcar(directory, nimages)
        data, trajectories = process_vasp_data(raw_data, directory, nimages, load_structs)
        return cls(data, trajectories)


def parse_outcar(directory, nimages):
    filename = directory + str(1).zfill(2) + "/OUTCAR"
    patterns = {"NTEMPER": r"parallel\stempering\sroutine\sentered\sNSTEP+=,\sNTEMPER+=\s+([\d\-.]+)\s+([\d\-.]+)",
                "Acceptance Ratio": r"Acceptance\sratio\sfor\sswaps" + "\s+([\d\-.]+)" * (nimages - 1),
                "old TEBEG": r"parallel\stempering\sold\sTEBEG" + "\s+([\d\-.]+)" * nimages,
                "new TEBEG": r"parallel\stempering\snew\sTEBEG" + "\s+([\d\-.]+)" * nimages,
                "old TOTEN": r"parallel\stempering\sold\sTOTEN" + "\s+([\d\-.]+)" * nimages,
                "attempt": r"attempting\sswapping\s+([\d\-.]+)",
                "dT": r"1/T1-1/T2\s+([\d\-.]+)",
                "dE": r"E1\s+-E2\s+([\d\-.]+)",
                "random_num": r"random\s+([\d\-.]+)",
                "success": r"swapping\s+([\d\-.]+)",
                "failure": r"noswapping\s+([\d\-.]+)"}

    return regrep(filename, patterns)


def process_vasp_data(raw_data, directory, nimages, load_structs=False):
    processing_data = deepcopy(raw_data)
    data = {}

    # Identify swap steps
    step_lines = [line[1] for line in processing_data["NTEMPER"]]
    swap_lines = [line for line_data, line in processing_data["Acceptance Ratio"]]
    data["nswap"] = [bisect.bisect_right(step_lines, i) for i in swap_lines]

    # clean up raw data
    for line in processing_data["attempt"]:
        swap_bool_lines = [line[1] for line in processing_data["success"]]
        index = bisect.bisect_left(swap_bool_lines, line[1])
        del processing_data["success"][index]

    success_fail = []
    failure_lines = [line[1] for line in processing_data["failure"]]
    for line in processing_data["success"]:
        if line[1] in failure_lines:
            success = False
        else:
            success = True
        index = bisect.bisect_right(swap_lines, line[1])
        success_fail.append((data["nswap"][index], line[0][0], success))

    swap_steps = {}
    for attempt in zip(success_fail, raw_data["dT"], raw_data["dE"], raw_data["random_num"]):
        step_num = attempt[0][0]
        if step_num not in swap_steps.keys():
            swap_steps[step_num] = []
        doc = {"dT": float(attempt[1][0][0]),
               "dE": float(attempt[2][0][0]),
               "random_number": float(attempt[3][0][0]),
               "swap_id": int(attempt[0][1]),
               "swap_bool": attempt[0][2]}
        swap_steps[step_num].append(doc)
    data["swap_steps"] = swap_steps

    temps = [[int(float(i)) for i in processing_data["old TEBEG"][0][0]]]
    for i in range(len(data["nswap"]) - 1):
        temps.extend(
            [[int(float(i)) for i in processing_data["new TEBEG"][i][0]]] * (data["nswap"][i + 1] - data["nswap"][i]))
    nsteps = 100
    temps.extend([[int(float(i)) for i in processing_data["new TEBEG"][i][0]]] * (nsteps - data["nswap"][i] - 2))
    image_swap_path = np.transpose(temps)
    data["image_trajectories"] = temps

    trajectories = {}
    if load_structs:
        structures = [None for i in range(nimages)]
        for i, temp in enumerate(temps[0]):
            filename = directory + str(i + 1).zfill(2) + "/XDATCAR"
            xdc = Xdatcar(filename)
            structures[i] = xdc.structures

        get_index = lambda x: temps[0].index(x)
        structures_reconstructed = [[] for i in range(nimages)]
        for image in range(nimages):
            for i, temp in enumerate(image_swap_path[image]):
                try:
                    structures_reconstructed[image].append(structures[get_index(temp)][i])
                except:
                    print(i)

        for temp, structures_image in zip(temps[0], structures):
            trajectories[temp] = Trajectory.from_structures(structures_image)

    # attempt_step = [bisect.bisect_right(step_lines, attempt[1]) for attempt in raw_data["attempt"]]
    #         print(attempt_step)
    #         for swap in data["success"][::2]:
    #             print(swap)
    data["old TOTEN"] = [i[0] for i in processing_data.get("old TOTEN", [])]
    data["acceptance"] = [i[0] for i in processing_data.get("Acceptance Ratio", [])]
    data["temperatures"] = temps[0]
    return data, trajectories
