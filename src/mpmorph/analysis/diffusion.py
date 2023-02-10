"""
This module contains new classes for obtaining
Diffusion and Activation Barrier calculations
from MD calculations.
"""

__author__ = "Muratahan Aykol <maykol@lbl.gov>"

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core import Element
from pymatgen.io.vasp import Xdatcar
from scipy import integrate, stats


class Diffusion(object):
    """
    Robust calculation of diffusion coefficients with different statistical
    analysis techniques:
        - Block averaging (default)
        - Jackknife (to be implemented)
        - Bootstrap (to be implemented)

    Args:
        structures: (list) list of Structures
        corr_t: (float) correlation time (in terms of # of steps).
            Each time origin will be this many steps apart.
        block_l: (int)  defines length of a block in terms of corr_t.
            (block_t = block_l * corr_t)
        t_step: (float) time-step in MD simulation. Defaults to 2.0 fs.
        l_lim: (int) this many time-steps are skipped in MSD while fitting D.
            i.e. approximate length of ballistic and cage regions. Defaults to 50.
        skip_first: (int) this many initial time-steps are skipped. Defaults to 0.
        ci: (float) confidence interval desired estimating the mean D of population.
    """

    def __init__(
        self, structures, corr_t, block_l, t_step=2.0, l_lim=50, skip_first=0, ci=0.95
    ):
        self.structures = structures
        self.abc = self.structures[0].lattice.abc
        self.natoms = len(self.structures[0])
        self.skip_first = skip_first
        self.total_t = len(self.structures)
        self.corr_t = corr_t
        self.l_lim = l_lim
        self.t_step = t_step
        self.block_l = block_l
        self.ci = ci
        self.msds = None
        self.vel_matrix = None
        self.vacfs = None
        self.scaling_factor = 0.1 / self.t_step  # conv. to cm2/s

    @property
    def n_origins(self):
        n = int((self.total_t - self.block_t - self.skip_first) / self.corr_t + 1)
        if n <= 0:
            raise ValueError("Too many blocks for the correlation time")
        return n

    @property
    def block_t(self):
        return self.block_l * self.corr_t

    def _getd(self, el):
        md = [np.zeros(self.structures[0].frac_coords.shape)]
        for i in range(self.skip_first + 1, self.total_t):
            dx = self.structures[i].frac_coords - self.structures[i - 1].frac_coords
            dx -= np.round(dx)
            md.append(dx)

        self.md = np.array(md) * self.abc

        # remove other elements from the rest of the calculations
        s = set(self.structures[0].indices_from_symbol(el))
        self.md = np.delete(
            self.md, [x for x in list(range(self.natoms)) if x not in s], 1
        )

        msds = []
        for i in range(self.n_origins):
            su = np.square(
                np.cumsum(
                    self.md[i * self.corr_t : i * self.corr_t + self.block_t], axis=0
                )
            )
            msds.append(np.mean(su, axis=1))
        self.msds = msds

    def plot_block_msds(self):
        for i in self.msds:
            plt.plot(i)

    def getD(self, el):
        """
        Method to calculate diffusion coefficient(s) of the given element (el).
        """

        if type(el) == type(Element("Li")):
            el = el.name

        self._getd(el)
        D = [[], [], []]
        for i in self.msds:
            for j in range(3):
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    np.arange(self.l_lim, self.block_t), i[:, j][self.l_lim :]
                )
                D[j].append(slope / 2.0)
        D = np.array(D) * self.scaling_factor
        self.D_blocks = D

        alpha = 1.0 - self.ci
        tn = stats.t.ppf(1.0 - alpha / 2.0, len(self.D_blocks) - 1) / np.sqrt(
            len(self.D_blocks)
        )

        if tn == "nan":
            tn = 1
        self.D_i = np.mean(D, axis=1)
        self.D_i_std = np.std(D, axis=1) * tn
        self.D_avg = np.sum(self.D_i) / 3.0
        self.D_avg_std = np.std(np.sum(D, axis=0) / 3.0) * tn
        return self.D_dict

    @property
    def D_dict(self):
        D_dict = {}
        dirs = ["Dx", "Dy", "Dz"]
        D_dict.update(dict(zip(dirs, self.D_i)))
        D_dict.update(dict(zip([s + "_std" for s in dirs], self.D_i_std)))
        D_dict.update({"D": self.D_avg, "D_std": self.D_avg_std})
        return D_dict

    @property
    def tao(self):
        tao_dict = {}
        for k, v in self.D_dict.items():
            if "_std" not in k:
                tao_dict[k] = 1.0 / v * self.scaling_factor
        return tao_dict

    def autocorrelation(self):
        "to be implemented"
        pass

    def get_v(self, el):
        # Make copy of structures
        _structures = [structure.copy() for structure in self.structures]

        # Find unneccessary elements and delete
        prune_els = []
        for specie in _structures[0].species:
            if specie != el:
                prune_els.append(specie)
        for structure in _structures:
            structure.remove_species(prune_els)
        _structures_sites = [structure.sites for structure in _structures]

        # Iterate through each site through each timestep and find velocity
        vel_matrix = [
            [0 for y in range(len(_structures) - 1)]
            for x in range(len(_structures[0].sites))
        ]
        for i in range(len(vel_matrix)):
            for j in range(len(vel_matrix[0])):
                vel_matrix[i][j] = (
                    _structures_sites[j][i].distance(_structures_sites[j + 1][i])
                    / self.t_step
                )
        self.vel_matrix = vel_matrix
        return

    def get_v_vector(self, el):
        # Make copy of structures
        _structures = [structure.copy() for structure in self.structures]

        # Find unneccessary elements and delete
        prune_els = []
        for specie in _structures[0].species:
            if specie != el:
                prune_els.append(specie)
        for structure in _structures:
            structure.remove_species(prune_els)
        _structures_sites = [structure.sites for structure in _structures]

        # Iterate through each site through each timestep and find velocity
        vel_matrix = [
            [[0, 0, 0] for y in range(len(_structures) - 1)]
            for x in range(len(_structures[0].sites))
        ]
        for i in range(len(vel_matrix)):
            for j in range(len(vel_matrix[0])):
                dist_x = _structures_sites[j][i].x - _structures_sites[j + 1][i].x
                if dist_x > _structures[i].lattice.a / 2:
                    dist_x = (_structures[i].lattice.a - np.abs(dist_x)) * (
                        -1 * np.sign(dist_x)
                    )
                dist_y = _structures_sites[j][i].y - _structures_sites[j + 1][i].y
                # if dist_y > _structures[i].lattice.b/2:
                #    dist_y = (_structures[i].lattice.b-np.abs(dist_y))*(-1*np.sign(dist_y))
                dist_z = _structures_sites[j][i].z - _structures_sites[j + 1][i].z
                # if dist_z > _structures[i].lattice.c/2:
                #    dist_z = (_structures[i].lattice.c-np.abs(dist_z))*(-1*np.sign(dist_z))

                vel_matrix[i][j][0] = dist_x / self.t_step
                vel_matrix[i][j][1] = dist_y / self.t_step
                vel_matrix[i][j][2] = dist_z / self.t_step

        self.vel_matrix = vel_matrix
        return

    def green_kubo_D(self, el):
        self.get_v(el)
        # Get velocity autocorrelation function for each site
        vacfs = []
        for site_vel in self.vel_matrix:
            _vacf = np.correlate(site_vel, site_vel, "full")
            vacfs.append(_vacf)
        self.vacfs = vacfs
        D = []
        for vacf in vacfs:
            D.append(integrate.simps(vacf))
        return D


class Activation(object):
    def __init__(self, D_t):
        self.D_t = D_t
        self.Q = None
        self.intercept = None
        self.Q_std = None

    def LS(self):
        self.x = np.array([1 / float(t[0]) for t in self.D_t])
        self.y = np.array([np.log(t[1]["D"]) for t in self.D_t])
        self.yerr = np.array(
            [
                [
                    -np.log((t[1]["D"] - t[1]["D_std"]) / t[1]["D"]),
                    np.log((t[1]["D"] + t[1]["D_std"]) / t[1]["D"]),
                ]
                for t in self.D_t
            ]
        )
        (
            self.Q,
            self.intercept,
            self.r_value,
            self.p_value,
            self.std_err,
        ) = stats.linregress(self.x, self.y)
        self.Q *= -1
        return self.Q

    def ODR(self):
        if not self.Q:
            self.LS()
        import scipy.odr

        def fit_func(p, t):
            return p[0] * t + p[1]

        Model = scipy.odr.Model(fit_func)
        Data = scipy.odr.RealData(self.x, self.y, sy=np.mean(self.yerr, axis=1))
        Odr = scipy.odr.ODR(Data, Model, [-self.Q, self.intercept])
        Odr.set_job(fit_type=2)
        self.output = Odr.run()
        self.Q, self.intercept = -self.output.beta[0], self.output.beta[1]
        self.Q_std = self.output.sd_beta[0]
        self.intercept_std = self.output.sd_beta[1]
        return self.Q, self.Q_std

    def plot(self, title=None, annotate=True, el="", **kwargs):
        # fig = plt.figure()

        line = np.polyval([-self.Q, self.intercept], self.x)
        tx = str(int(np.rint(self.Q)))
        if self.Q_std:
            tx += "$\pm${}".format(str(int(np.rint(self.Q_std))))
        c = kwargs.get("color", "")
        plt.plot(
            self.x * 1000,
            line,
            c + "-",
        )
        plt.errorbar(
            self.x * 1000,
            self.y,
            yerr=self.yerr.T,
            label="Q[{}]: ".format(el) + tx + " K",
            **kwargs,
        )
        plt.ylabel("ln(D cm$^2$/s)", fontsize=15)
        plt.xlabel("1000/T K$^{-1}$", fontsize=15)

        if annotate:
            plt.annotate(
                "Q: " + tx + " K",
                xy=(0.98, 0.95),
                xycoords="axes fraction",
                fontsize=14,
                horizontalalignment="right",
                verticalalignment="top",
            )
        if title:
            plt.title = title
            # return fig

    @classmethod
    def from_run_paths(
        cls, p, T, el, corr_t, block_l, t_step=2.0, l_lim=50, skip_first=0
    ):
        D_t = []
        for t in range(len(p)):
            xdatcar = Xdatcar(p[t])
            d = Diffusion(
                xdatcar.structures,
                corr_t=corr_t,
                block_l=block_l,
                t_step=t_step,
                l_lim=l_lim,
                skip_first=skip_first,
            )
            D_t.append([T[t], d.getD(el)])
        return cls(D_t)
