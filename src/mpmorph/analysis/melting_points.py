import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import mean_squared_error
import numpy as np
import matplotlib.pyplot as plt
import math
from abc import ABC, abstractmethod
import numpy.polynomial.polynomial as poly
from numpy.polynomial import Polynomial as P

class AbstractMeltingPointEstimator(ABC):

    def plot(self, ts, vs, plot_title=None):
        fig, axs = self._plot_ts_vs(ts, vs)
        Tm = self.estimate(ts, vs)
        axs.plot([Tm, Tm], [min(vs), max(vs)], color="r")

        if plot_title is None:
            axs.set_title("Volume vs Temperature by Polynomial Fit")
        else:
            axs.set_title(plot_title)
        
        return fig, axs

    @abstractmethod
    def estimate(self, temps, vols):
        pass

    def _plot_ts_vs(self, ts, vs):
        fig, axs = plt.subplots()
        axs.scatter(ts, vs)
        axs.set_xlabel("Temperature (K)")
        axs.set_ylabel("Volume (A^3)")
        return fig, axs        

class MeltingPointEnsembleEstimator(AbstractMeltingPointEstimator):

    def __init__(self, estimators):
        self._estimators = estimators
    
    def estimate(self, temps, vols):
        tm_estimates = [e.estimate(temps, vols) for e in self._estimators]
        return np.mean(tm_estimates)

    def plot(self, ts, vs, plot_title = None):
        _, axs = self._plot_ts_vs(ts, vs)
        tm = self.estimate(ts, vs)
        axs.plot([tm, tm], [np.min(vs), np.max(vs)])
    
    def plot_all_estimates(self, ts, vs):
        _, axs = self._plot_ts_vs(ts, vs)
        for e in self._estimators:
            tm = e.estimate(ts, vs)
            axs.plot([tm, tm], [np.min(vs), np.max(vs)], label=e.name)

        avg = self.estimate(ts, vs)
        axs.plot([avg , avg ], [np.min(vs), np.max(vs)], label="Mean Estimate")
        axs.set_title("Ensemble Tm Estimates")
        axs.legend()



class MeltingPointClusterEstimator(AbstractMeltingPointEstimator):

    name: str = "Clustering"

    def plot_clusters(self, ts, vs, plot_title=None):
        points = np.array(list(zip(ts, vs)))
        cluster1, cluster2 = self._get_clusters(points)
        fig, axs = plt.subplots()
        axs.scatter(*cluster1)
        axs.scatter(*cluster2)
        axs.set_xlabel("Temperature (K)")
        axs.set_ylabel("Volume (A^3)")
        Tm = self.estimate(ts, vs)
        axs.plot([Tm, Tm], [min(vs), max(vs)], color="r")

        if plot_title is None:
            axs.title("Volume vs Temperature by Clustering")
        else:
            axs.title(plot_title)
        
        return fig, axs

    def estimate(self, temps, vols):
        points = np.array(list(zip(temps, vols)))
        cluster1, cluster2 = self._get_clusters(points)
        if min(cluster1[0]) < min(cluster2[0]):
            solid_range = cluster1[0]
            liquid_range = cluster2[0]
        else:
            solid_range = cluster2[0]
            liquid_range = cluster1[0]

        return np.mean([max(solid_range), min(liquid_range)])

    def _get_clusters(self, points):
        clustering = AgglomerativeClustering(n_clusters=2).fit(points)
        cluster1 = points[np.argwhere(clustering.labels_ == 1).squeeze()].T
        cluster2 = points[np.argwhere(clustering.labels_ == 0).squeeze()].T
        return cluster1, cluster2

class MeltingPointTrisectionEstimator(AbstractMeltingPointEstimator):

    def _get_fit_error_total(self, xs, ys):
        slope, intercept, r_value, p_value, std_err = linregress(xs, ys)
        y_pred = intercept + slope * np.array(xs)
        err = np.sum(np.abs(y_pred - ys))        
        return slope, intercept, err

    def _unzip_pts(self, pts):
        ts = [pt[0] for pt in pts]
        vs = [pt[1] for pt in pts]
        return ts, vs

    def _plot_pts(self, pts):
        ts, vs = self._unzip_pts(pts)
        plt.scatter(ts,vs, color='grey')


    def _split_pts(self, pts, x1, x2):
        set1 = [pt for pt in pts if pt[0] < x1]
        set2 = [pt for pt in pts if pt[0] > x1 and pt[0] < x2]
        set3 = [pt for pt in pts if pt[0] > x2]
        return set1, set2, set3

    def _plot_split(self, pts, x1, x2):
        set1, set2, set3 = self._plit_pts(pts, x1, x2)
        self._plot_pts(set1)
        self._plot_pts(set2)
        self._plot_pts(set3)

    def _get_linear_ys(self, m, b, xs):
        return [m * x + b for x in xs]

    def _plot_fit_line(self, m, b, xs):
        fit_ys = self._get_linear_ys(m, b, xs)
        plt.plot(xs, fit_ys)

    def _plot_fits(self, pts, x1, x2):
        set1, set2, set3 = self._split_pts(pts, x1, x2)
        for dset in [set1, set2, set3]:
            xs, ys = self._unzip_pts(dset)
            m, b, err = self._get_fit_error_total(xs, ys)
            self._plot_pts(dset)
            self._plot_fit_line(m, b, xs)


    def _find_best_trisection(self, points, min_x = None, max_x = None, min_window_size = 100, step_size = 50):
        lowest_err = math.inf
        
        xs, ys = self._unzip_pts(points)
        if min_x is None:
            min_x = math.floor(np.min(xs))

        if max_x is None:
            max_x = math.ceil(np.max(xs))
        
        for pt1 in range(min_x + min_window_size, max_x - 2 * min_window_size, step_size):
            for pt2 in range(pt1 + min_window_size, max_x - min_window_size, step_size):
                set1, set2, set3 = self._split_pts(points, pt1, pt2)
                errs_total = 0
                errs = []
                try:
                    for dset in [set1, set2, set3]:
                        xs, ys = self._unzip_pts(dset)
                        m, b, err = self._get_fit_error_total(xs, ys)      
                        errs.append(err)
                        errs_total += err ** 2

                    errs_total = math.sqrt(errs_total)

                    if errs_total < lowest_err:
                        lowest_err = errs_total
                        best_pt1 = pt1
                        best_pt2 = pt2
                except:
                    print("problem encountered")

        return best_pt1, best_pt2

    def _estimate_melting_pt(self, points):
        pt1, pt2 = self._find_best_trisection(points, step_size=100)
        # print(f'Coarse guess: {pt1, pt2}')
        pt1, pt2 = self._find_best_trisection(points, pt1 - 150, pt2 + 150, step_size = 5)
        # print(f'Fine guess: {pt1, pt2}')
        return pt1, pt2, (pt2 - pt1) / 2 + pt1
    
    def estimate(self, temps, vols):
        pts = list(zip(temps, vols))
        pt1, pt2, tm = self._estimate_melting_pt(pts)
        return tm
    
    def plot(self, ts, vs, plot_title=None):
        pts = list(zip(ts, vs))
        pt1, pt2, tm = self._estimate_melting_pt(pts)
        self._plot_fits(pts, pt1, pt2)


class MeltingPointBisectionEstimator(AbstractMeltingPointEstimator):

    def plot_best_split(self, temps, vols):
        split_idx = self.get_best_split(temps, vols)
        fig, axs = self._plot_split(temps, vols, split_idx)
        Tm = self.estimate(temps, vols)
        axs.plot([Tm, Tm], [min(vols), max(vols)], color="r")
        return fig, axs

    def estimate(self, temps, vols):
        best_split_idx = self.get_best_split(temps, vols)
        return np.mean([temps[best_split_idx], temps[best_split_idx - 1]])

    def _split_dset(self, pts, split_idx):
        return pts[0:split_idx], pts[split_idx:]

    def _assess_splits(self, xs, ys):
        dset_size = len(xs)
        buffer = max(round(dset_size / 10), 3)

        pt_idxs = list(range(buffer + 1, len(xs) - buffer - 1))
        errs = []
        for idx in pt_idxs:
            _, _, _, _, total_err = self._get_split_fit(xs, ys, idx)
            errs.append(total_err)

        return list(zip(pt_idxs, errs))

    def _get_linear_ys(self, m, b, xs):
        return [m * x + b for x in xs]

    def _plot_split(self, xs, ys, split_idx):
        m1, b1, m2, b2, _ = self._get_split_fit(xs, ys, split_idx)
        leftxs, rightxs = self._split_dset(xs, split_idx)
        leftys, rightys = self._split_dset(ys, split_idx)

        fig, axs = plt.subplots()

        left_fit_ys = self._get_linear_ys(m1, b1, leftxs)
        axs.scatter(leftxs, leftys)
        axs.plot(leftxs, left_fit_ys)


        right_fit_ys = self._get_linear_ys(m2, b2, rightxs)
        axs.scatter(rightxs, rightys)
        axs.plot(rightxs, right_fit_ys)

        axs.title("Volume vs Temperature (w/ best fits by Slope Method")
        axs.xlabel("Temperature (K)")
        axs.ylabel("Equil. Volume (cubic Angstroms)")
        return fig, axs

    def get_best_split(self, xs, ys):
        split_errs = self._assess_splits(xs, ys)
        errs = [pt[1] for pt in split_errs]
        idxs = [pt[0] for pt in split_errs]
        best_split_idx = idxs[np.argmin(errs)]
        return best_split_idx


class MeltingPointSlopeRMSEEstimator(MeltingPointBisectionEstimator):

    name: str = "RMSE Bisection"

    def _get_fit_error(self, xs, ys):
        slope, intercept, r_value, p_value, std_err = linregress(xs, ys)
        y_pred = intercept + slope * np.array(xs)
        err = mean_squared_error(y_true=ys, y_pred=y_pred, squared=False)        
        return slope, intercept, err

    def _get_split_fit(self, xs, ys, split_idx):
        leftx, rightx = self._split_dset(xs, split_idx)
        lefty, righty = self._split_dset(ys, split_idx)

        lslope, lintercept, r_value, p_value, std_err = linregress(leftx, lefty)
        left_y_pred = lintercept + lslope * np.array(leftx)
        lefterr = mean_squared_error(y_true=lefty, y_pred=left_y_pred, squared=False)

        rslope, rintercept, r_value, p_value, std_err = linregress(rightx, righty)
        right_y_pred = rintercept + rslope * np.array(rightx)
        righterr = mean_squared_error(y_true=righty, y_pred=right_y_pred, squared=False)

        combined_err = math.sqrt(lefterr**2 + righterr**2)
        combined_err = lefterr + righterr
        return lslope, lintercept, rslope, rintercept, combined_err


class MeltingPointSlopeStdErrEstimator(MeltingPointBisectionEstimator):

    name: str = "StdErr Bisection"

    def _get_split_fit(self, xs, ys, split_idx):
        leftx, rightx = self._split_dset(xs, split_idx)
        lefty, righty = self._split_dset(ys, split_idx)

        leftfit = linregress(leftx, lefty)
        lefterr = leftfit.stderr

        rightfit = linregress(rightx, righty)
        righterr = rightfit.stderr

        combined_err = math.sqrt(lefterr**2 + righterr**2)
        combined_err = lefterr + righterr
        return (
            leftfit.slope,
            leftfit.intercept,
            rightfit.slope,
            rightfit.intercept,
            combined_err,
        )

class MeltingPointPolyFitEstimator(AbstractMeltingPointEstimator):

    name: str = "Polynomial Fit"

    def plot_fit(self, ts, vs, plot_title=None):
        fig, axs = super().plot(ts, vs, plot_title=plot_title)
        p = self._polyfit(ts, vs)
        fit_ts = p(vs)
        axs.plot(fit_ts, vs, color='r')
        return fig, axs

    def estimate(self, temps, vols):
        p = self._polyfit(temps, vols)
        second = p.deriv(2)
        melting_v = second.roots()[0]
        melting_t = p(melting_v)
        return melting_t
    
    def _polyfit(self, temps, vols):
        coefs = poly.polyfit(vols, temps, 3)
        p = P(coefs)
        return p