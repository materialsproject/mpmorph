from copy import deepcopy
from multiprocessing import Pool, RawArray
import numpy as np
from scipy.optimize import curve_fit
import bisect
import os
from monty.json import MSONable
from monty.functools import lru_cache
from monty.re import regrep
from pymatgen import Structure
from pymatgen.io.vasp.outputs import Xdatcar


class Trajectory(MSONable):
    def __init__(self, base_structure, disp, lattices, time_step=2):
        self.time_step = time_step
        self.base_structure = base_structure
        self.base_frac_coords = base_structure.frac_coords
        self.disp = np.array(disp)
        self.lattices = lattices
        self.index = 0
        self.structure = base_structure

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
            self.lattices = np.concatenate((self.lattices, trajectory.lattices))

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

        return cls(structure, disp, l)

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
             }
        return d

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, d["displacements"], d["lattices"])


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
                plt.plot([i, i], [min(self.get_temperatures()), max(self.get_temperatures())],
                         '--k', linewidth=1)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("time step (#)", fontsize=18)
        plt.ylabel("Temperature", fontsize=18)
        plt.show()

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "data": self.data,
             "trajectories": [trajectory.as_dict()
                              for trajectory in self.trajectories]
             }
        return d

    @classmethod
    def from_dict(cls, d):
        trajectories = [Trajectory.from_dict(trajectory_dict)
                        for trajectory_dict in d["trajectories"]]
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


# A global dictionary storing the variables passed from the initializer.
var_dict = {}


def init_worker(shared_d, d_shape):
    var_dict['shared_d'] = shared_d
    var_dict['d_shape'] = d_shape


class MSD(MSONable):
    def __init__(self, t, msd, std_dev, temperature, time_step, sites, site_multiplicity):
        self.t = t
        self.msd = msd
        self.std_dev = std_dev
        self.temperature = temperature
        self.time_step = time_step
        self.sites = sites
        self.site_multiplicity = site_multiplicity

    def average_species(self):
        avg_msd = self.msd.copy()
        avg_msd_sites = self.sites.copy()
        avg_std = self.std_dev.copy()
        site_mult = self.site_multiplicity.copy()
        # Ensure the msd is not already averaged over species:
        if len(self.sites) < np.shape(avg_msd)[0]:
            _avg_msd = []
            _avg_std = []
            _avg_sites = {}
            _site_mult = {}
            for i, (key, items) in enumerate(avg_msd_sites.items()):
                _avg_msd.append(np.mean(np.array(avg_msd)[items], axis=0))
                #                 _avg_std.append(np.std(np.array(avg_msd)[items], axis=0))

                # Compute as sum of normally distributed random variables
                # (may not be accurate given that x, y, z are not random):
                _avg_std.append(np.sqrt(np.divide(np.sum(np.square(np.array(avg_std)[items]), axis=0), len(items))))
                _avg_sites[key] = [i]
                _site_mult[key] = [len(items)]
            avg_msd = _avg_msd
            avg_std = _avg_std
            avg_msd_sites = _avg_sites
            site_mult = _site_mult
        return MSD(self.t, avg_msd, avg_std, temperature=self.temperature, time_step=self.time_step,
                   sites=avg_msd_sites, site_multiplicity=site_mult)

    def avg_msds(self, msds):
        """
        Average over multiple msds

        :param msds:
        :return:
        """
        # Ensure the msd is already averaged over species:
        all_msds = [self.average_species()]
        for msd in msds:
            all_msds.append(msd.average_species())
        # if len(self.sites) < np.shape(self.msd)[0]:
        #             self.average(['specie', 'vector'])
        avg_msds = []
        avg_stds = []
        site_mult = {}
        for key in self.sites.keys():
            # TODO: Need to ensure all trajectories are the same length, else truncate to the shortest

            # CalculateAverage
            msd_list = [msd.msd[msd.sites[key][0]] * msd.site_multiplicity[key]
                        for msd in all_msds]
            msd_array = np.stack(msd_list)
            msd_sum = np.sum(msd_array, axis=0)
            total_sites = np.sum([msd.site_multiplicity[key] for msd in all_msds])
            msd_avg = msd_sum / total_sites

            # Calculate std deviation
            std_list = [msd.std_dev[msd.sites[key][0]] for msd in all_msds]
            std_array = np.stack(std_list)
            std_dev = np.sqrt(np.divide(np.sum(np.square(std_array), axis=0), len(std_list)))

            site_mult[key] = total_sites
            avg_msds.append(msd_avg)
            avg_stds.append(std_dev)
        avg_msd_sites = self.sites
        return MSD(self.t, avg_msds, avg_stds, temperature=self.temperature,
                   time_step=self.time_step, sites=avg_msd_sites,
                   site_multiplicity=site_mult)

    @classmethod
    def get_msd_from_structures(cls, structures, temperature, specie, x_vals=None, skip=1,
                                time_step=2):
        traj = Trajectory.from_structures(structures)
        return cls.get_msd_from_displacements(traj.structure, traj.disp, temperature,
                                              specie, x_vals, skip, time_step)

    @classmethod
    def get_msd_from_displacements(cls, structure, displacements, temperature, specie,
                                   x_vals=None, skip=1, time_step=2):
        """

        Args:
            structure:
            displacements: Numpy array with shape  [site, time step, axis]
            temperature:
            x_vals:
            skip:
            time_step:

        Returns:

        """
        framework_indices = list(structure.indices_from_symbol(specie))
        framework_disp = displacements[framework_indices]
        drift = np.average(framework_disp, axis=0)[None, :, :]

        # drift corrected position and reshape to [time step, site, axis]
        dc = displacements - drift
        dc = np.swapaxes(dc, 0, 1)[::skip]
        nsteps, nions, dim = dc.shape
        if type(x_vals) not in [np.ndarray, list]:
            # Default to evenly spaced points in log-space (Makes for prettier log-log plots)
            run_length = nsteps
            x_vals = np.unique(
                np.logspace(0, np.log10(run_length - 1),
                            int(np.log10(run_length - 1) * 100), base=10, dtype=int))
        # Calculate MSD at each x_val using multiprocessing
        shared_d = RawArray('d', nsteps * nions * dim)
        shared_d_np = np.frombuffer(shared_d, dtype=np.float64).reshape((nsteps, nions, dim))
        np.copyto(shared_d_np, dc)
        # from multiprocessing import Pool
        with Pool(processes=32, initializer=init_worker,
                  initargs=(shared_d, (nsteps, nions, dim))) as pool:
            results = pool.map(process_chunk, [(h, i) for h, i in enumerate(x_vals)])
        # Extract data from results
        #         msd = np.zeros((len(x_vals), np.shape(displacements)[1], 3))
        #         msd_std = np.zeros((len(x_vals), np.shape(displacements)[1], 3))
        msd = np.zeros((len(x_vals), nions))
        msd_std = np.zeros((len(x_vals), nions))
        for result in results:
            index = result[0]
            msd[index] = result[1]
            msd_std[index] = result[2]
        # Swap the axes so the array is structured as [atom][time_step]
        msd = np.swapaxes(msd, 0, 1)
        msd_std = np.swapaxes(msd_std, 0, 1)
        t = np.multiply(x_vals, skip * time_step)
        # Get the sites and site multiplicity
        species = structure.composition.elements
        sites = dict([(str(el), []) for el in species])
        site_mult = dict([(str(el), []) for el in species])
        for i, el in enumerate(structure.species):
            sites[str(el)].append(i)
        for key in sites.keys():
            sites[key] = np.array(sites[key])
            site_mult[key] = np.ones(len(sites[key]))
        return cls(t, msd, msd_std, time_step=time_step, sites=sites,
                   site_multiplicity=site_mult, temperature=temperature)

    def plot(self, axes='None', atomic_index=0, **kwargs):
        from matplotlib import pyplot as plt
        plt.figure(figsize=kwargs.get('figsize', [7, 5]), dpi=kwargs.get('dpi', 150))

        ax = plt.axes()
        if axes in ['semilogx', 'log']:
            ax.set_xscale("log")
        if axes in ['semilogy', 'log']:
            ax.set_yscale("log")
        for key in self.sites:
            ax.errorbar(self.t, self.msd[self.sites[key][atomic_index]],
                        yerr=self.std_dev[self.sites[key][atomic_index]], label=f'{key}')
        ax.set_ylabel(r'MSD ($\AA{}^2$)', fontsize=18)
        ax.set_xlabel('time (ps)', fontsize=18)

        plt.legend(fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        return plt


class DiffusionAnalyzer(MSONable):
    def __init__(self, msds):
        """
        :param msds: (list) of MSD objects
        """
        self.msds = msds
        self.msds_dict = {}
        for msd in msds:
            try:
                self.msds_dict[msd.temperature].append(msd)
            except:
                self.msds_dict[msd.temperature] = [msd]
        self.msds_avg = {}
        for key, msd_list in self.msds_dict.items():
            self.msds_avg[key] = msd_list[0].avg_msds(msd_list[1:])
        self.temps = sorted(self.msds_avg.keys(), reverse=True)
        self.elements = sorted(list(self.msds_avg[self.temps[0]].sites.keys()), key = lambda x: self.msds_avg[self.temps[0]].sites[x])
        d = []
        d_std = []
        for temp in self.temps:
            # TODO: Check if the MSD meets a minimum of 1 (or maybe a higher cutoff)
            _d, _d_std = self.fit_msd(self.msds_avg[temp])
            d.append(_d)
            d_std.append(_d_std)
        self.d = np.array(d)
        self.d_std = np.array(d_std)
        self.fit_data = {}
        self.room_temp_d = {}
        self.room_temp_d_std = {}
        self.plot_data = {}
        for i, el in enumerate(self.elements):
            x = self.temps
            y = self.d[:, i]
            yerr = self.d_std[:, i]
            popt, pcov = curve_fit(arrhenius_eqn, x, y, sigma=yerr, absolute_sigma=True)
            self.fit_data[el]=(popt, pcov)
            # Calculate room temperature value
            room_t = 298.15
            d_fit = arrhenius_eqn(room_t, *popt)
            # Calculate error on extrapolated value
            sigma_a = np.sqrt(np.diag(pcov))[0]
            sigma_e = np.sqrt(np.diag(pcov))[1]
            sigma_d = np.sqrt((sigma_a ** 2 * (d_fit / popt[0]) ** 2) + sigma_e ** 2 * (d_fit / extrap_t) ** 2)
            #set
            self.room_temp_d[el] = d_fit
            self.room_temp_d_std[el] = sigma_d
#     print(fr'{d_fit:.2E} +- {sigma_d:.2E} cm2/s')

    def get_plot_data(self):
        plot_data = {}
        for i, el in enumerate(self.elements):
            popt, pcov = self.fit_data[el]
            x = self.temps
            x_scaled = [1000 / t for t in x]
            y = self.d[:, i]
            fit_line_x = x.copy()
            fit_line_x.append(298.15)
            fit_line_x_scaled = [1000 / t for t in fit_line_x]
            fit_line_y = arrhenius_eqn(fit_line_x, *popt)
            plot_data[el] = (x_scaled, y, fit_line_x_scaled, fit_line_y)
        return plot_data

    def plot(self):
        from matplotlib import pyplot as plt
        plot_data = self.get_plot_data()
        plt.figure(figsize=[7, 5], dpi=150)
        ax = plt.axes()
        ax.set_yscale("log")
        colors = ['blue', 'green', 'orange']
        for i, el in enumerate(self.elements):
            x_scaled, y, fit_line_x_scaled, fit_line_y = plot_data[el]
            ax.errorbar(x_scaled, y, self.d_std[:, i], color=colors[i], label=str(el))
            ax.errorbar(fit_line_x_scaled, fit_line_y, fmt='--', color=colors[i], label=None)
            ax.errorbar([1000/298.15], [self.room_temp_d[el]], [self.room_temp_d_std[el]], fmt='x', label=None)
        plt.legend()
        ax.set_ylabel(r'D ($cm^2/s$)', fontsize=18)
        ax.set_xlabel('1000/T ($K^{-1}$)', fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        return plt

    def fit_msd(self, msd):
        d_arr = []
        d_std_arr = []
        min_t = np.where(msd.t > 1)[0][0] - 1
        max_t = np.where(msd.t > msd.t[-1] * .5)[0][0] - 1
        for i, el in enumerate(self.elements):
        # for i in range(len(msd.msd)):
            popt, pcov = curve_fit(linear_eqn, msd.t[min_t:max_t], msd.msd[i][min_t:max_t],
                                   sigma=msd.std_dev[i][min_t:max_t], absolute_sigma=True)
            d_arr.append(popt[0] * 1E-4)
            d_std_arr.append(np.sqrt(np.diag(pcov))[0] * 1E-4)
        return d_arr, d_std_arr


def arrhenius_eqn(t, a, b):
    return a * np.exp(-b / t)


def linear_eqn(x, m, b):
    return m * x + b


def process_chunk(data):
    """
    This function computes the time-averaged point for one time

    :param data: a tuple of (int, int) corresponding to index and dt for this chunk
    :return:
    """

    h = data[0]
    i = data[1]
    a = np.frombuffer(var_dict['shared_d']).reshape(var_dict['d_shape'])
    b = np.subtract(a[i:], a[:-i])
    c = np.square(b)
    del a, b
    msd = np.mean(np.sum(c, axis=-1), axis=0)
    msd_std = np.std(np.sum(c, axis=-1), axis=0)
    del c
    return h, msd, msd_std