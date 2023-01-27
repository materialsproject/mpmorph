import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import mean_squared_error
import numpy as np
import matplotlib.pyplot as plt
import math

class MeltingPointClusterAnalyzer():

    def _get_clusters(self, points):
        clustering = AgglomerativeClustering(n_clusters=2).fit(points)
        cluster1 = points[np.argwhere(clustering.labels_ == 1).squeeze()].T
        cluster2 = points[np.argwhere(clustering.labels_ == 0).squeeze()].T
        return cluster1, cluster2

    def plot_vol_vs_temp(self, ts, vs, plot_title = None):
        points = np.array(list(zip(ts, vs)))
        cluster1, cluster2 = self._get_clusters(points)
        plt.scatter(*cluster1)
        plt.scatter(*cluster2)
        plt.xlabel("Temperature (K)")
        plt.ylabel("Volume (A^3)")
        Tm = self.estimate_melting_temp(ts, vs)
        plt.plot([Tm, Tm], [min(vs), max(vs)], color='r')
      
        if plot_title is None:
            plt.title("Volume vs Temperature by Clustering")
        else:
            plt.title(plot_title)
    
    def estimate_melting_temp(self, temps, vols):
        points = np.array(list(zip(temps, vols)))
        cluster1, cluster2 = self._get_clusters(points)
        if min(cluster1[0]) < min(cluster2[0]):
            solid_range = cluster1[0]
            liquid_range = cluster2[0]
        else:
            solid_range = cluster2[0]
            liquid_range = cluster1[0]

        return np.mean([max(solid_range), min(liquid_range)])

class MeltingPointSlopeAnalyzer():

    def split_dset(self, pts, split_idx):
        return pts[0:split_idx], pts[split_idx:]

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
        plt.title("Volume vs Temperature (w/ best fits by Slope Method")
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
        Tm = self.estimate_melting_temp(temps, vols)
        print(Tm)
        plt.plot([Tm, Tm], [min(vols), max(vols)], color='r')        


    def estimate_melting_temp(self, temps, vols):
        best_split_idx = self.get_best_split(temps, vols)
        return np.mean([temps[best_split_idx], temps[best_split_idx - 1]])

class MeltingPointSlopeRMSEAnalyzer(MeltingPointSlopeAnalyzer):

    def get_split_fit(self, xs, ys, split_idx):
        leftx, rightx = self.split_dset(xs, split_idx)
        lefty, righty = self.split_dset(ys, split_idx)
        
        lslope, lintercept, r_value, p_value, std_err = linregress(leftx, lefty)
        left_y_pred = lintercept + lslope * np.array(leftx)
        lefterr = mean_squared_error(y_true=lefty, y_pred=left_y_pred, squared=False)

        rslope, rintercept, r_value, p_value, std_err = linregress(rightx, righty)
        right_y_pred = rintercept + rslope * np.array(rightx)
        righterr = mean_squared_error(y_true=righty, y_pred=right_y_pred, squared=False)
        
        combined_err = math.sqrt(lefterr ** 2 + righterr ** 2)
        combined_err = lefterr + righterr
        return lslope, lintercept, rslope, rintercept, combined_err

class MeltingPointSlopeStdErrAnalyzer(MeltingPointSlopeAnalyzer):

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
