"""
This module contains new classes for obtaining
Diffusion and Activation Barrier calculations
from MD calculations.
"""

__author__ = 'Muratahan Aykol <maykol@lbl.gov>'

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pymatgen.io.vasp import Xdatcar


class Diffusion(object):
    """
    Robust calculation of diffusion coefficients with different statistical analysis techniques:
        - Block averaging (default)
        - Jackknife (to be implemented)
        - Bootstrap (to be implemented)

    Args:
        structures: (list) list of Structures
        corr_t: (float) correlation time. Each time origin will be this many steps apart.
        block_l: (int)  defines length of a block in terms of corr_t. (block_t = block_l * corr_t)
        t_step: (float) time-step in MD simulation. Defaults to 2.0 fs.
        l_lim: (int) this many time-steps are skipped in MSD while fitting D. I.e. approximate length of
            ballistic and cage regions. Defaults to 50.
        skip_first: (int) this many initial time-steps are skipped. Defaults to 0.
    """

    def __init__(self, structures, corr_t, block_l, t_step=2.0, l_lim=50, skip_first=0):
        self.structures = structures
        self.abc = self.structures[0].lattice.abc
        self.natoms = len(self.structures[0])
        self.skip_first = skip_first
        self.total_t = len(self.structures)
        self.corr_t = corr_t
        self.l_lim = l_lim
        self.t_step = t_step
        self.block_l = block_l  #
        self.msds = None
        self.scaling_factor = 0.1 / self.t_step  # conv. to cm2/s

    @property
    def n_origins(self):
        n = (self.total_t - self.block_t - self.skip_first) / self.corr_t + 1
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
        self.md = np.delete(self.md, [x for x in list(range(self.natoms)) if x not in s], 1)

        msds = []
        for i in range(self.n_origins):
            su = np.square(np.cumsum(self.md[i * self.corr_t: i * self.corr_t + self.block_t], axis=0))
            msds.append(np.mean(su, axis=1))
        self.msds = msds

    def plot_block_msds(self):
        for i in self.msds:
            plt.plot(i)

    def getD(self, el):
        """
        Method to calculate diffusion coefficient(s) of the given element (el).
        """
        self._getd(el)
        D = [[], [], []]
        for i in self.msds:
            for j in range(3):
                slope, intercept, r_value, p_value, std_err = \
                    stats.linregress(np.arange(self.l_lim, self.block_t), i[:, j][self.l_lim:])
                D[j].append(slope / 2.0)
        D = np.array(D) * self.scaling_factor
        self.D_blocks = D
        self.D_i = np.mean(D, axis=1)
        self.D_i_std = np.std(D, axis=1)
        self.D_avg = np.sum(self.D_i) / 3.0
        self.D_avg_std = np.std(np.sum(D, axis=0) / 3.0)
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


class Activation(object):
    def __init__(self, D_t):
        self.D_t = D_t
        self.Q = None
        self.intercept = None
        self.Q_std = None

    def LS(self):
        self.x = np.array([1 / float(t[0]) for t in self.D_t])
        self.y = np.array([np.log(t[1]["D"]) for t in self.D_t])
        self.yerr = np.array([[-np.log((t[1]["D"] - t[1]["D_std"]) / t[1]["D"]),
                               np.log((t[1]["D"] + t[1]["D_std"]) / t[1]["D"])
                               ] for t in self.D_t])
        self.Q, self.intercept, self.r_value, self.p_value, self.std_err = \
            stats.linregress(self.x, self.y)
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
        return self.Q, self.Q_std

    def plot(self, title=None):
        fig = plt.figure()
        plt.errorbar(self.x * 1000, self.y, yerr=self.yerr.T, fmt='o')
        line = np.polyval([-self.Q, self.intercept], self.x)
        plt.plot(self.x * 1000, line, 'k-')
        plt.ylabel("ln(D cm$^2$/s)", fontsize=15)
        plt.xlabel("1000/T K$^{-1}$", fontsize=15)
        tx = str(int(np.rint(self.Q)))
        if self.Q_std:
            tx += "$\pm${}".format(str(int(np.rint(self.Q_std))))
        plt.annotate("Q: " + tx + " K", xy=(0.98, 0.95), xycoords='axes fraction', fontsize=14,
                     horizontalalignment='right', verticalalignment='top')
        if title:
            plt.title = title
        return fig

    @classmethod
    def from_run_paths(cls, p, T, el, corr_t, block_l, t_step=2.0, l_lim=50, skip_first=0):
        D_t = []
        for t in range(len(p)):
            xdatcar = Xdatcar(p[t])
            d = Diffusion(xdatcar.structures, corr_t=corr_t, block_l=block_l,
                          t_step=t_step, l_lim=l_lim, skip_first=skip_first)
            D_t.append([T[t], d.getD(el)])
        return cls(D_t)