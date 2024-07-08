"""
Developed by Dr. Angshuman Deb with Degenkolb Engineers.
It provides tools for analyzing hysteresis curves and determine energy dissipation ratio.
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from shapely.geometry import LineString, Polygon
from scipy.signal import find_peaks
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mplPolygon
from scipy.stats import uniform
import piecewise_regression
from matplotlib import colors


def display_ticks(ax):
    ax.minorticks_on()
    ax.grid(True, which="major", alpha=0.6)
    ax.grid(True, which="minor", alpha=0.6)
    return ax


def sync_axes(axs, **kwargs):
    x = kwargs.get("x", True)
    y = kwargs.get("y", True)
    z = kwargs.get("z", True)

    """
    Synchronize the x and y-axis limits across a list of axes.

    :param axs: List of matplotlib axes objects.
    """
    xmin, xmax, ymin, ymax = float('inf'), float('-inf'), float('inf'), float('-inf')
    zmin = float('inf')
    zmax = float('-inf')

    # Find the global min and max limits
    for ax in axs:
        xlimits = ax.get_xlim()
        ylimits = ax.get_ylim()

        try:
            zlimits = ax.get_zlim()
        except AttributeError:
            zlimits = None

        xmin = min(xmin, xlimits[0])
        xmax = max(xmax, xlimits[1])
        ymin = min(ymin, ylimits[0])
        ymax = max(ymax, ylimits[1])

        if zlimits is not None:
            zmin = min(zmin, zlimits[0])
            zmax = max(zmax, zlimits[1])

    # Set global limits to all axes
    for ax in axs:
        if x:
            ax.set_xlim([xmin, xmax])
        if y:
            ax.set_ylim([ymin, ymax])
        try:
            if z:
                ax.set_zlim([zmin, zmax])
        except AttributeError:
            pass


def create_uniform_dist(number, scale_factor=0.2):
    """
    - scale_factor (float, optional): The proportion of the number to determine the magnitude of variability. Defaults to 0.2.
    - if number is 0, then the distribution will be built around (-1, 1)
    """
    magnitude = abs(number) * scale_factor
    if magnitude == 0:
        magnitude = 1
    lower_bound = number - magnitude
    upper_bound = number + magnitude
    dist = uniform(loc=lower_bound, scale=upper_bound - lower_bound)
    return dist


def piecewise_linear(x, *params):
    # Break up params into two lists: breakpoints and slopes
    b = params[:(len(params) - 1) // 2]
    k = params[(len(params) - 1) // 2:]

    # Define the function for each piece
    def piece_func(x_piece, seg_num):
        return k[seg_num] * x_piece + sum([(k[itr_sn] - k[itr_sn + 1]) * b[itr_sn] for itr_sn in range(seg_num)])

    condlist = [x < b[0]] + [(x >= b[itr_b]) & (x < b[itr_b + 1]) for itr_b in range(len(b) - 1)] + [x >= b[-1]]
    # Apply the correct function for each x
    y = np.zeros_like(x)
    for itr_line_seg in range(len(k)):
        y[condlist[itr_line_seg]] = piece_func(x[condlist[itr_line_seg]], itr_line_seg)
    return y


def resample_array(original_array, m):
    """

    :param original_array:
    :param m: number of points in resampled array
    :return:
    """
    n = len(original_array)
    # Check if m is smaller than n
    if m > n:
        raise ValueError("m should be smaller than n.")
    spacing = (n - 1) / (m - 1)
    fill_indices = [int(round(i * spacing)) for i in range(m)]
    resampled_array = [original_array[i] for i in fill_indices]

    return np.array(resampled_array)


class HysteresisCurve:
    def __init__(self, d, f, n_seg_env, **kwargs):
        """
        :param d: displacement history (positive is assumed direction of loading, and a cycle will be unloading, then loading)
        :param f: force history
        :param n_seg_env: number of piecewise linear segments in positive f-d envelope
        :param kwargs:
        """
        self.labels = kwargs.get("labels", ('D', 'F'))
        self.prominence = kwargs.get("prominence", None)
        self.d = np.array(d)
        self.f = np.array(f)
        self.n_seg_env = n_seg_env
        self.env_params = None

        d_m = []
        f_m = []
        max_d = 0.
        for i, _d in enumerate(self.d):
            if _d >= max_d:
                d_m.append(_d)
                f_m.append(self.f[i])
                max_d = _d
        self.d_m = np.array(d_m)
        self.f_m = np.array(f_m)

    def find_peaks_and_valleys(self, prominence=None):
        if prominence is None:
            prominence = self.prominence
        peaks, _ = find_peaks(self.d, prominence=prominence)
        valleys, _ = find_peaks(-self.d, prominence=prominence)
        valleys = valleys[valleys > peaks[0]]
        return peaks, valleys

    def see_peaks_and_valleys(self, prominence):
        peaks, valleys = self.find_peaks_and_valleys(prominence)
        fig, axs = plt.subplots(2, 1, sharex="all", sharey="all", constrained_layout=True)
        steps = np.arange(len(self.d))
        axs[0].plot(steps, self.d, 'b-')
        axs[1].plot(steps, self.d, 'b-')
        axs[0].plot(steps[peaks], self.d[peaks], 'r^')
        axs[1].plot(steps[valleys], self.d[valleys], 'rv')
        for ax in axs:
            display_ticks(ax)
        fig.supxlabel("Step number")
        fig.supylabel(self.labels[0])
        axs[0].set_title("Peaks")
        axs[1].set_title("Valleys")

    def resample(self, n):
        """

        :param n: number of points in each unloading/reloading branch (Total number of points in 1 cycle = 2 * n)
        :return:
        """
        ubs = self.unloading_branches()
        rbs = self.reloading_branches()
        if len(ubs) > 0:
            d = [*resample_array(self.d[:ubs[0][0]], n)]
            f = [*resample_array(self.f[:ubs[0][0]], n)]
        else:
            d = [*resample_array(self.d, n)]
            f = [*resample_array(self.f, n)]
        for i in range(len(ubs)):
            d.extend(resample_array(self.d[ubs[i][0]:ubs[i][1]], n))
            d.extend(resample_array(self.d[rbs[i][0]:rbs[i][1]], n))
            f.extend(resample_array(self.f[ubs[i][0]:ubs[i][1]], n))
            f.extend(resample_array(self.f[rbs[i][0]:rbs[i][1]], n))
        return HysteresisCurve(np.array(d), np.array(f), self.n_seg_env)

    def get_init_slope_estimate(self):
        return self.get_env_params()[self.n_seg_env - 1]

    def unloading_branches(self):
        peaks, valleys = self.find_peaks_and_valleys()
        return list(zip(peaks, valleys))

    def reloading_branches(self):
        peaks, valleys = self.find_peaks_and_valleys()
        peaks = peaks[1:]
        if len(peaks) == len(valleys) - 1:
            ub = self.unloading_branches()
            if len(ub) > 0 and ub[-1][-1] == valleys[-1]:
                peaks = np.hstack((peaks, len(self.d)))
        return list(zip(valleys, peaks))

    def hysteresis_polygons(self):
        u_branches = self.unloading_branches()
        r_branches = self.reloading_branches()
        cycles = len(u_branches)
        m1 = self.get_init_slope_estimate()
        env_params = self.get_env_params()
        yield_pt = env_params[0]
        polygons = []
        for cycle in range(cycles):
            u_branch = u_branches[cycle]
            r_branch = r_branches[cycle]
            f_u = self.f[u_branch[0]:u_branch[1] + 1].copy()
            d_u = self.d[u_branch[0]:u_branch[1] + 1].copy()
            f_r = self.f[r_branch[0]:r_branch[1] + 1].copy()
            d_r = self.d[r_branch[0]:r_branch[1] + 1].copy()
            if abs(d_u[0]) < yield_pt and abs(d_u[-1]) < yield_pt:
                polygons.append((Polygon(), Polygon()))
            else:
                polygons.append(HysteresisCurve._hysteresis_polygons(d_u, f_u, d_r, f_r, m1))
        return polygons

    def edr(self):
        polygons = self.hysteresis_polygons()
        edr = []
        for hyst_poly, epp_poly in polygons:
            if hyst_poly.area == 0 or epp_poly.area == 0:
                edr.append(0.)
            # elif hyst_poly.area / epp_poly.area > 1:
            #     edr.append(0.)
            else:
                edr.append(hyst_poly.area / epp_poly.area)
        return np.array(edr)

    def plot(self, plot_cycle_num=None, **kwargs):
        axs = kwargs.get('axs', None)
        color = kwargs.get('color', 'b')
        plot_env = kwargs.get('plot_env', False)

        alpha = 0.7
        if isinstance(color, str):
            color = colors.to_rgba(color)[:-1]

        if axs is None:
            _, ax0 = plt.subplots(1, 1, constrained_layout=True)
            _, ax1 = plt.subplots(1, 1, constrained_layout=True)
            if plot_cycle_num is not None:
                _, ax2 = plt.subplots(1, 1, constrained_layout=True)
                axs = [ax0, ax1, ax2]
            else:
                axs = [ax0, ax1]
        axs[0].plot(self.d, self.f, '-', color=color + (alpha,))

        edr = self.edr()
        cycles = np.arange(0, len(edr))
        if len(cycles) > 0:
            ml, sl, bl = axs[1].stem(cycles, edr)
            ml.set_markeredgecolor("none")
            ml.set_markerfacecolor(color + (alpha,))
            sl.set_color(color + (alpha,))
            bl.set_color(color + (alpha,))
            axs[1].set_ylim(-0.1, np.maximum(1., 1.1 * np.max(edr)))

        if plot_cycle_num is not None:
            if np.isscalar(plot_cycle_num):
                plot_cycle_num_list = [plot_cycle_num]
            else:
                plot_cycle_num_list = plot_cycle_num
            for plot_cycle_num in plot_cycle_num_list:
                try:
                    polygons = self.hysteresis_polygons()
                    hyst_polygon = mplPolygon(np.array(polygons[plot_cycle_num][0].exterior.coords), closed=True, facecolor='red',
                                              edgecolor='none', alpha=0.25)
                    epp_polygon = mplPolygon(np.array(polygons[plot_cycle_num][1].exterior.coords), closed=True, facecolor='none',
                                             edgecolor='k', alpha=0.5)
                    axs[0].add_patch(hyst_polygon)
                    axs[0].add_patch(epp_polygon)
                    ml, sl, bl = axs[1].stem(cycles[plot_cycle_num], [edr[plot_cycle_num]])
                    ml.set_markeredgecolor("none")
                    ml.set_markerfacecolor(color + (1.,))
                    ml.set_marker('*')
                    ml.set_markersize(20)
                    sl.set_color(color + (1.,))
                    bl.set_color(color + (1.,))

                    hyst_polygon = mplPolygon(np.array(polygons[plot_cycle_num][0].exterior.coords), closed=True,
                                              facecolor='red',
                                              edgecolor=color, alpha=0.25)
                    epp_polygon = mplPolygon(np.array(polygons[plot_cycle_num][1].exterior.coords), closed=True,
                                             facecolor='none',
                                             edgecolor='k', alpha=0.5)
                    axs[2].add_patch(hyst_polygon)
                    axs[2].add_patch(epp_polygon)
                except IndexError:
                    continue

        if plot_env:
            d_e = np.linspace(min(self.d_m), max(self.d_m), 100)
            f_e = piecewise_linear(d_e, *self.get_env_params())
            axs[0].plot(d_e, f_e, 'k--')
            axs[0].plot(-d_e, -f_e, 'k--')
            if plot_cycle_num:
                axs[2].plot(d_e, f_e, 'k--')
                axs[2].plot(-d_e, -f_e, 'k--')

        display_ticks(axs[0])
        display_ticks(axs[1])
        if plot_cycle_num:
            display_ticks(axs[2])

        axs[0].set_xlabel(self.labels[0])
        axs[0].set_ylabel(self.labels[1])

        axs[1].set_xlabel('cycle')
        axs[1].set_ylabel('EDR')

        if plot_cycle_num:
            axs[2].set_xlabel(self.labels[0])
            axs[2].set_ylabel(self.labels[1])

        abs_xlim = max(abs(np.array(axs[0].get_xlim())))
        abs_ylim = max(abs(np.array(axs[0].get_ylim())))
        axs[0].set_xlim((-abs_xlim, abs_xlim))
        axs[0].set_ylim((-abs_ylim, abs_ylim))
        if plot_cycle_num:
            sync_axes([axs[0], axs[2]])
        return axs

    def get_env_params(self):
        """
        Recommended for monotonic force-deformation
        """
        if self.env_params is None:
            x_data, y_data = self.d_m, self.f_m

            def fun(p):
                return np.sum((y_data - piecewise_linear(x_data, *p)) ** 2)

            n_breakpoints = self.n_seg_env - 1
            prf = piecewise_regression.Fit(list(x_data), list(y_data), n_breakpoints=n_breakpoints)
            est = prf.get_results()['estimates']
            bps = []
            for i in range(1, n_breakpoints + 1):
                bps.append(est[f'breakpoint{i}']['estimate'])

            inds = []

            for i in range(len(bps) + 1):
                if i == 0:
                    inds.append(x_data < bps[i])
                elif i == len(bps):
                    inds.append(x_data > bps[i - 1])
                else:
                    inds.append((x_data > bps[i - 1]) & (x_data < bps[i]))

            ks = []
            for ind in inds:
                ks.append(np.polyfit(x_data[ind], y_data[ind], 1)[0])

            dists = [create_uniform_dist(item) for item in bps + ks]
            p0s = np.array([dist.rvs(10, random_state=311) for dist in dists]).T

            res_list = []
            p_list = []
            for p0 in p0s:
                temp = minimize(
                    fun, x0=p0,
                    bounds=[
                        u.support() for u in dists
                    ],
                )
                res_list.append(temp.fun)
                p_list.append(temp.x)
            self.env_params = p_list[np.argmin(res_list)]
        return self.env_params

    @staticmethod
    def find_intersection(x1, y1, m1, x2, y2, m2):
        x = (y2 - m2 * x2 - y1 + m1 * x1) / (m1 - m2)
        y = y1 + m1 * (x - x1)
        return x, y

    @staticmethod
    def find_other_vertices(x1, y1, x3, y3, m1, m2):
        x2, y2 = HysteresisCurve.find_intersection(x1, y1, m2, x3, y3, m1)
        x4, y4 = HysteresisCurve.find_intersection(x1, y1, m1, x3, y3, m2)
        return (x2, y2), (x4, y4)

    @staticmethod
    def _hysteresis_polygons(d_u, f_u, d_r, f_r, m1):
        ind = d_r < d_u[0]
        d_r = d_r[ind]
        f_r = f_r[ind]
        ucl = LineString(zip(d_u, f_u))
        rcl = LineString(zip(d_r, f_r))
        intersection = rcl.intersection(ucl)
        if intersection.geom_type == 'Point':
            d_r = np.concatenate((d_r, np.array([d_u[0]])))
            f_r = np.concatenate((f_r, np.array([f_u[0]])))
            ucl = LineString(zip(d_u, f_u))
            rcl = LineString(zip(d_r, f_r))
            intersection = rcl.intersection(ucl)

        x_points = sorted(list(intersection.geoms), key=lambda p: (p.coords[0][0], p.coords[0][1]))
        if len(x_points) > 2:
            x_points = [x_points[0], x_points[-1]]
        other_points = HysteresisCurve.find_other_vertices(*x_points[0].coords[0], *x_points[1].coords[0], m1, 0.)
        epp_poly = Polygon(
            (x_points[0].coords[0], other_points[0], x_points[1].coords[0], other_points[1], x_points[0].coords[0]))

        idx_u = [np.argmin(np.linalg.norm(np.column_stack((d_u, f_u)) - point.coords[0], axis=1)) for point in x_points]
        idx_r = [np.argmin(np.linalg.norm(np.column_stack((d_r, f_r)) - point.coords[0], axis=1)) for point in x_points]

        for itr in range(len(x_points)):
            d_u[idx_u[itr]] = x_points[itr].coords[0][0]
            f_u[idx_u[itr]] = x_points[itr].coords[0][1]

            d_r[idx_r[itr]] = x_points[itr].coords[0][0]
            f_r[idx_r[itr]] = x_points[itr].coords[0][1]

        d_u = d_u[min(idx_u):max(idx_u) + 1]
        f_u = f_u[min(idx_u):max(idx_u) + 1]

        d_r = d_r[min(idx_r):max(idx_r) + 1]
        f_r = f_r[min(idx_r):max(idx_r) + 1]

        seg1 = np.column_stack((d_u, f_u))
        seg2 = np.column_stack((d_r[1:], f_r[1:]))
        try:
            hyst_poly = Polygon(np.concatenate((seg1, seg2)))
        except ValueError:
            hyst_poly = Polygon()
        return hyst_poly, epp_poly


if __name__ == "__main__":
    d_vals = np.array(
        [-2.06044411e-19, 4.00000000e-01, 8.00000000e-01, 1.10000000e+00,
         1.50000000e+00, 1.90000000e+00, 2.30000000e+00, 2.60000000e+00,
         3.00000000e+00, 3.40000000e+00, 3.50000000e+00, 2.70000000e+00,
         2.00000000e+00, 1.20000000e+00, 4.00000000e-01, -3.00000000e-01,
         -1.10000000e+00, -1.90000000e+00, -2.60000000e+00, -3.40000000e+00,
         -3.50000000e+00, -2.50000000e+00, -1.60000000e+00, -6.00000000e-01,
         3.00000000e-01, 1.30000000e+00, 2.20000000e+00, 3.20000000e+00,
         4.10000000e+00, 5.10000000e+00, 5.20000000e+00, 4.10000000e+00,
         2.90000000e+00, 1.80000000e+00, 6.00000000e-01, -5.00000000e-01,
         -1.70000000e+00, -2.80000000e+00, -4.00000000e+00, -5.10000000e+00,
         -5.20000000e+00, -3.70000000e+00, -2.10000000e+00, -6.00000000e-01,
         9.00000000e-01, 2.50000000e+00, 4.00000000e+00, 5.50000000e+00,
         7.10000000e+00, 8.60000000e+00, 8.70000000e+00, 6.80000000e+00,
         4.90000000e+00, 2.90000000e+00, 1.00000000e+00, -9.00000000e-01,
         -2.80000000e+00, -4.80000000e+00, -6.70000000e+00, -8.60000000e+00,
         -8.70000000e+00, -6.60000000e+00, -4.50000000e+00, -2.40000000e+00,
         -3.00000000e-01, 1.90000000e+00, 4.00000000e+00, 6.10000000e+00,
         8.20000000e+00, 1.03000000e+01, 1.04000000e+01, 8.10000000e+00,
         5.80000000e+00, 3.50000000e+00, 1.20000000e+00, -1.20000000e+00,
         -3.50000000e+00, -5.80000000e+00, -8.10000000e+00, -1.04000000e+01,
         -1.05000000e+01, -8.00000000e+00, -5.50000000e+00, -3.00000000e+00,
         -5.00000000e-01, 2.10000000e+00, 4.60000000e+00, 7.10000000e+00,
         9.60000000e+00, 1.21000000e+01]
    )
    f_vals = np.array(
        [-8.06023761e-26, 7.34016040e+00, 1.46796495e+01, 2.01838141e+01,
         2.75211834e+01, 3.48365937e+01, 4.20911136e+01, 4.73902542e+01,
         5.35760564e+01, 6.03991104e+01, 6.19185495e+01, 4.79748571e+01,
         3.57664780e+01, 2.17476856e+01, 7.30456754e+00, -5.53148738e+00,
         -2.02042711e+01, -3.48163260e+01, -4.65012937e+01, -6.00741136e+01,
         -6.13104865e+01, -4.39486731e+01, -2.83257056e+01, -1.08526393e+01,
         5.53183145e+00, 2.36480288e+01, 3.93493017e+01, 5.66728212e+01,
         6.54560002e+01, 6.83328835e+01, 6.86175240e+01, 4.96368539e+01,
         2.95371906e+01, 1.28446928e+01, 2.64539993e-01, -9.25866403e+00,
         -3.05788302e+01, -4.95801496e+01, -6.45254034e+01, -6.82810439e+01,
         -6.85776856e+01, -4.27907994e+01, -1.70619410e+01, -2.97987445e-01,
         1.61778538e+01, 4.22452429e+01, 6.04677517e+01, 6.82926943e+01,
         7.34231968e+01, 7.77556704e+01, 7.80390556e+01, 4.56775826e+01,
         1.80925376e+01, 4.57327231e+00, -2.07827393e+00, -1.65971197e+01,
         -4.73209913e+01, -6.53336154e+01, -7.21745217e+01, -7.77308730e+01,
         -7.80151142e+01, -4.23836595e+01, -1.40666470e+01, -2.64533132e+00,
         4.18405911e+00, 3.33460851e+01, 5.94795111e+01, 6.95625550e+01,
         7.62908050e+01, 8.01459865e+01, 8.01128926e+01, 4.05636988e+01,
         1.30949412e+01, 2.92408129e+00, -3.92376404e+00, -2.64374889e+01,
         -5.40820538e+01, -6.79709923e+01, -7.54984225e+01, -8.01058058e+01,
         -8.00727120e+01, -3.62727916e+01, -1.02268830e+01, -6.68479335e-01,
         6.40894677e+00, 3.20048337e+01, 5.80810806e+01, 7.00278994e+01,
         7.77646937e+01, 7.95243929e+01]
    )
    # fd = np.loadtxt("InputForceDisp.txt", delimiter=',')

    curve = HysteresisCurve(d_vals, f_vals, n_seg_env=2)

    # first use hc.see_peaks_and_valleys(prominence=p) for different values of p to determine the adequate prominence
    curve.see_peaks_and_valleys(prominence=0.55)
    # then set
    curve.prominence = 0.55
    curve.plot()
    curve.plot([2, 3], plot_env=False)
