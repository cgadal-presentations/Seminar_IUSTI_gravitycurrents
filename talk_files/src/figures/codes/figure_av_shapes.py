import os
import sys

import cmocean as cmo
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic
from uncertainties import ufloat
from uncertainties import unumpy as unp


def _complex_shape(Z, a, b):
    prefact = 1 / (2 * np.pi * (a**2 - 2 * a * b + 4 * b**2))
    term1 = 3 * a * (2 - a) * np.log(1 + 2 * b * Z / a)
    term2 = 3 * (2 * b**2 - a * b + a) * np.log(1 - Z) + (2 * a**2 - a * b + 2 * b**2 - 3 * a) * np.log(1 - Z**3)
    term3 = 2 * np.sqrt(3) * (2 * b**2 + a * b - 4 * b + a) * (np.arctan((2 * Z + 1) / np.sqrt(3)) - np.pi / 6)
    return prefact * (term1 + term2 + term3)


def _real_shape(R, a, b, Hc, Hs):
    Z = R * np.exp(1j * np.pi / 3)
    C = _complex_shape(Z, a=np.sqrt(3) * (2 / np.pi) ** (1 / 3), b=1.4316)
    x = np.real(C)
    y = np.imag(C)
    return Hc * x, -Hc * y + Hs


def Benjamin_shape(Hc, Hs, a=np.sqrt(3) * (2 / np.pi) ** (1 / 3), b=1.4316):
    R = np.linspace(0.0001, 100, 1000)
    x, y = _real_shape(R, a, b, Hc, Hs)
    interpolated_function = interp1d(x, y, kind="cubic", fill_value="extrapolate")
    return interpolated_function


def Benjamin_shape_fit(x, Hc, Hs):
    a = np.sqrt(3) * (2 / np.pi) ** (1 / 3)
    b = 1.4316
    xt, yt = _real_shape(100, a, b, Hc, Hs)
    return np.where(x < xt, Benjamin_shape(Hc, Hs)(x), yt)


def compute_av_shape(ds, ind_tmin, ind_tmax):
    x = np.concatenate([lx - np.nanmax(lx) for lx in ds.variables["contour time series (x)"][ind_tmin:ind_tmax]])
    y = np.concatenate([ly - np.nanmin(ly) for ly in ds.variables["contour time series (y)"][ind_tmin:ind_tmax]])
    #
    all_points = np.vstack([x, y]).T
    mask = np.isnan(all_points[:, 0]) | np.isnan(all_points[:, 1])
    all_points = np.array([all_points[:, 0][~mask], all_points[:, 1][~mask]]).T
    #
    bins = [
        int(np.abs(all_points[:, 0]).max() * 1e2 / ds["px_per_cm"][:].data),
        int(np.abs(all_points[:, 1]).max() * 1e2 / ds["px_per_cm"][:].data),
    ]
    # h, xedges, yedges = np.histogram2d(all_points[:, 0], all_points[:, 1], bins=bins)
    # xcenters = (xedges[1:] + xedges[:-1]) / 2
    # ycenters = (yedges[1:] + yedges[:-1]) / 2
    #
    av_shape, bin_edges, binnumber = binned_statistic(
        all_points[:, 0], all_points[:, 1], statistic="mean", bins=bins[0]
    )
    av_shape_std, bin_edges, binnumber = binned_statistic(
        all_points[:, 0], all_points[:, 1], statistic="std", bins=bins[0]
    )
    av_shape_counts, bin_edges, binnumber = binned_statistic(
        all_points[:, 0], all_points[:, 1], statistic="count", bins=bins[0]
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return av_shape, av_shape_std, av_shape_counts, bin_centers


# path data
path_gen = "/media/cyril/LaCie_orang/petit_canal"
round_manip = "round_spring2022"
results_dir = "Processing/Results"
ref = "sand80m_H19/"  # directoy to process
path_total = os.path.join(path_gen, round_manip, results_dir, ref)
par_dir = os.path.join(path_total, "../../scripts/image_processing/par_files")
data_dir = os.path.join(path_total, "shape")

# Loading data
Position_processed = np.load(
    os.path.join(path_total, "nose_position", "Position_processed.npy"),
    allow_pickle=True,
).item()
# Loading shapes
SHAPES = np.load(os.path.join(data_dir, "av_shapes/Av_shapes.npy"), allow_pickle=True).item()
SHAPES_props = np.load(os.path.join(data_dir, "av_shapes_log/Shape_logs_props.npy"), allow_pickle=True).item()

# Loading initial parameters
L_reservoir = ufloat(9.9, 0.2)  # [cm]
W_reservoir = ufloat(19.4, 0.2)  # [cm]
#
Parameters = np.load(os.path.join(path_total, "Initial_parameters.npy"), allow_pickle=True).item()

runs = sorted(Position_processed.keys())
Volume_fraction = np.array([Parameters[run]["Volume_fraction"] for run in runs])
runs_sorted, phi_sorted = np.array([[run, phi] for phi, run in sorted(zip(Volume_fraction, runs))]).T
ind = -1
runs_sorted, phi_sorted = runs_sorted[1:], phi_sorted[1:]
#

cmap = cmo.cm.haline_r
log_phi = np.log10(unp.nominal_values(phi_sorted))
colors = cmap((log_phi - log_phi.min()) / (0.8 * log_phi.max() - log_phi.min()))
colors[:, -1] = 1

examples = ["run01", "run03", "run11", "run19", "run15", "run02"]

for ifig in range(3):
    print("draw figures")
    figsize = (1.3 * quarto.regular_fig_width, 0.75 * quarto.regular_fig_height)
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize)

    for run, color, phi in zip(runs_sorted, colors, phi_sorted):
        H0 = Parameters[run]["V_reservoir"] / W_reservoir / L_reservoir  # [cm]
        if run in examples:
            ax.plot(
                SHAPES[run]["xcenters"] / L_reservoir.n,
                unp.nominal_values(SHAPES[run]["shape"]) / H0.n,
                label=f"${100 * phi.n:.1f}$",
                color=color,
            )
    if ifig > 0:
        ax.axvline(-1, color="w", ls="--")
    fig.legend(title=r"$\phi~[\%]$", loc="outside right center")
    #
    if ifig > 1:
        xplot = np.linspace(-12, 0, 100)
        (a,) = ax.plot(
            xplot,
            Benjamin_shape_fit(
                np.abs(xplot),
                0.9,
                0,
            ),
            color="w",
            lw=2,
            label="Benjamin's shape",
        )
        leg1 = ax.legend(handles=[a], loc="lower center")

    ax.set_xlabel(r"Distance behind front, $(x - x_{\rm f})/l_{0}$")
    ax.set_ylabel(r"Height, $h/h_{0}$")
    ax.set_xlim(left=-8)
    ax.set_ylim(bottom=0)

    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
