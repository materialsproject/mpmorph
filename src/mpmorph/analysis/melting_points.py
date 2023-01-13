import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

import math

class MeltingPointAnalyzer():

    def split_dset(self, pts, split_idx):
        return pts[0:split_idx], pts[split_idx:]

    def get_split_fit(self, xs, ys, split_idx):
        leftx, rightx = self.split_dset(xs, split_idx)
        lefty, righty = self.split_dset(ys, split_idx)
        
        leftfit = linregress(leftx, lefty)
        lefterr = leftfit.stderr
        
        rightfit = linregress(rightx, righty)
        righterr = rightfit.stderr
        
        combined_err = math.sqrt(lefterr ** 2 + righterr ** 2)
        combined_err = lefterr + righterr
        return leftfit.slope, leftfit.intercept, rightfit.slope, rightfit.intercept, combined_err


    def assess_splits(self, xs, ys):
        dset_size = len(xs)
        buffer = max(round(dset_size / 10), 3)

        pt_idxs = list(range(buffer + 1, len(xs) - buffer - 1))
        errs = []
        for idx in pt_idxs:
            _, _, _, _, total_err = self.get_split_fit(xs, ys, idx)
            errs.append(total_err)
        
        return list(zip(pt_idxs, errs))

    def get_linear_ys(self, m, b, xs):
        return [m * x + b for x in xs]

    def plot_split(self, xs, ys, split_idx):
        m1, b1, m2, b2, _ = self.get_split_fit(xs, ys, split_idx)
        leftxs, rightxs = self.split_dset(xs, split_idx)
        leftys, rightys = self.split_dset(ys, split_idx)

        left_fit_ys = self.get_linear_ys(m1, b1, leftxs)

        plt.scatter(leftxs, leftys)
        plt.plot(leftxs, left_fit_ys)
        plt.title("Volume vs Temperature (w/ best fits)")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Equil. Volume (cubic Angstroms)")

        right_fit_ys = self.get_linear_ys(m2, b2, rightxs)

        plt.scatter(rightxs, rightys)
        plt.plot(rightxs, right_fit_ys)
        
    def get_best_split(self, xs, ys):
        split_errs = self.assess_splits(xs, ys)
        errs = [pt[1] for pt in split_errs]
        idxs = [pt[0] for pt in split_errs]
        best_split_idx = idxs[np.argmin(errs)]
        return best_split_idx
    
    def plot_vol_vs_temp(self, temps, vols):
        split_idx = self.get_best_split(temps, vols)
        self.plot_split(temps, vols, split_idx)

    def estimate_melting_temp(self, temps, vols):
        best_split_idx = self.get_best_split(temps, vols)
        return temps[best_split_idx]