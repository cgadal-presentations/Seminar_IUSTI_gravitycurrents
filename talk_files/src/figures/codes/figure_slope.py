import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
import template as tp
from matplotlib.colors import to_rgba
from matplotlib.patches import Ellipse
from models import Froude
from netCDF4 import Dataset

plt.rcParams["xtick.top"] = False


def ang2sin(ang):
    return np.sin(np.radians(ang))


def sin2ang(sin):
    return np.degrees(np.arcsin(sin))


# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, Fr, St, H0 = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
            d.variables["Fr"][:].data,
            d.variables["St"][:].data,
            d.variables["H0"][:].data,
        ]
        for d in datasets
    ]
).T

Ha = np.array(
    [d.variables["H_a"][:].data if "H_a" in d.variables.keys() else d.variables["H0"][:].data for d in datasets]
)
authors, particles = np.array(
    [
        [
            d.author,
            d.particle_type,
        ]
        for d in datasets
    ]
).T

# %% graphic specifications
# %% masks for plot
mask_phi = phi < tp.phi_c

# %% graphic vector for plots
alphas = np.ones_like(Fr)
alphas[~mask_phi] = 0.4

dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == "SedFoam"] = "s"

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
facecolors = np.array([to_rgba(c, a) for c, a in zip(facecolors, alphas)])

mask_nosuspended = (authors == "Marie Rastello") & (H0 / Ha < 1)
edgecolors = np.array([to_rgba("dimgray", a) for a in alphas])
edgecolors[mask_nosuspended] = np.array([to_rgba("tab:red", 0.4) for a in alphas[mask_nosuspended]])

zorders = np.vectorize(lambda dataset: tp.dataset_zorder2[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% theory
phi_m = 0.585
alpha_plot = np.linspace(0, 50, 500)
eta = 0
Fr_th = Froude(alpha_plot, eta, 7e4)

# %% figure

for ifig in range(3):
    figsize = (1.3 * quarto.regular_fig_width, quarto.regular_fig_height)
    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

    tp.mscatter(
        ang2sin(alpha)[plot_idxs],
        Fr[plot_idxs],
        ax=ax,
        m=markers[plot_idxs],
        facecolors=facecolors[plot_idxs],
        edgecolors=edgecolors[plot_idxs],
        lw=0.5,
    )

    # (line_birman,) = ax.plot(
    #     ang2sin(alpha_plot),
    #     Birman(alpha_plot),
    #     ls="--",
    #     # color="k",
    #     zorder=-8,
    #     lw=1,
    #     label="simulations of\nBirman et al. 2007",
    # )

    if ifig > 0:
        (line_theory,) = ax.plot(ang2sin(alpha_plot), Fr_th, ls="-", zorder=-8, lw=1, label="eq. (15)", color="w")

    if ifig > 1:
        ellipse = Ellipse((0, 0.35), 0.05, 0.55, facecolor="none", edgecolor="w")
        ax.add_patch(ellipse)

    ax.set_xlabel(r"$\sin \alpha$")
    ax.set_xlim(left=-0.012)

    secax = ax.secondary_xaxis("top", functions=(sin2ang, ang2sin))
    secax.set_xlabel(r"slope, $\alpha$ [deg.]")

    ax.set_ylabel(r"Slumping velocity, $\mathcal{F}r$")
    ax.set_ylim(0, 1.4)
    ax.set_xlim(right=ang2sin(48))

    leg1 = fig.legend(
        handles=tp.legend_datasets,
        ncol=1,
        title="Datasets",
        loc="outside right upper",
    )
    leg2 = fig.legend(handles=tp.legend_particles, ncol=1, title="Particles", loc="outside right center")

    fig.align_labels()

    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
