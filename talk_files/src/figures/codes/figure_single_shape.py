import glob
import os
import sys

import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from scipy.stats import binned_statistic

# path data
path_data = "/media/cyril/LaCie_orang/petit_canal/Dataset_shapes/Dataset/runs"

list_nc_files = glob.glob(os.path.join(path_data, "*.nc"))
g = 9.81  # [m/s2]

# nc_file = list_nc_files[50]
nc_file = list_nc_files[45]

ds = Dataset(nc_file)
h0 = ds.variables["H0"][:].data  # lock characteristic height, [m]
l0 = ds.variables["L0"][:].data  # lock characteristic width, [m]
rho_p = ds.variables["rho_p"][:].data  # particle velocity, [kg/m3]
rho_f = ds.variables["rho_f"][:].data  # lock fluid density, [kg/m3]
rho_a = ds.variables["rho_a"][:].data  # ambient fluid density, [kg/m3]
alpha = ds.variables["alpha"][:].data  # bottom slope, [deg.]
diam = ds.variables["d"][:].data  # grain size, [m]
phi = ds.variables["phi"][:].data
#
rho_c = rho_f + phi * (rho_p - rho_f)  # average lock density, [kg/m3]
gprime = g * (rho_c - rho_a) / rho_a  # specific gravity
u0 = np.sqrt(gprime * h0 * np.cos(np.radians(alpha)))
t_ad = l0 / u0
#
tmin, tmax = 3 * t_ad, 18 * t_ad
ind_tmin = np.argmin(np.abs(ds.variables["t"][:].data - tmin))
ind_tmax = np.argmin(np.abs(ds.variables["t"][:].data - tmax))

X, Y = [], []
for lx, ly in zip(
    ds.variables["contour time series (x)"][ind_tmin:ind_tmax],
    ds.variables["contour time series (y)"][ind_tmin:ind_tmax],
):
    X.append((lx - np.nanmax(lx)) / l0)
    Y.append((ly - np.nanmin(ly)) / h0)

# av shape
all_points = np.vstack([np.concatenate(X), np.concatenate(Y)]).T
mask = np.isnan(all_points[:, 0]) | np.isnan(all_points[:, 1])
all_points = np.array([all_points[:, 0][~mask], all_points[:, 1][~mask]]).T
#
print(np.abs(all_points[:, 0]).max())
bins = [
    int(np.abs(all_points[:, 0]).max() * 1e2 / ds["px_per_cm"][:].data),
    int(np.abs(all_points[:, 1]).max() * 1e2 / ds["px_per_cm"][:].data),
]
av_shape, bin_edges, binnumber = binned_statistic(all_points[:, 0], all_points[:, 1], statistic="mean", bins=bins[0])
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# figure
for ifig in range(2):
    figsize = (1.2 * quarto.regular_fig_width, 0.75 * quarto.regular_fig_height)
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize, sharex=True)
    for x, y in zip(X, Y):
        ax.plot(x, y, color="C1", alpha=0.075)
    if ifig > 0:
        ax.plot(bin_centers, av_shape, color="tab:blue", lw=2, label="average")
        ax.legend()
    ax.set_xlabel(r"Distance behind front, $(x - x_{\rm f})/l_{0}$")
    ax.set_ylabel(r"Height, $h/h_{0}$")
    ax.set_ylim(bottom=0)

    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
