import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
import template as tp
from models import Froude, Krieger_viscosity
from netCDF4 import Dataset

plt.rcParams["figure.constrained_layout.hspace"] = 0
plt.rcParams["figure.constrained_layout.h_pad"] = 0.005


# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, St, H0, L0, Fr = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
            d.variables["St"][:].data,
            d.variables["H0"][:].data,
            d.variables["L0"][:].data,
            d.variables["Fr"][:].data,
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
dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == "SedFoam"] = "s"

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
edgecolors = np.full_like(facecolors, "dimgray")

mask_nosuspended = (authors == "Marie Rastello") & (H0 / Ha < 1)
edgecolors[mask_nosuspended] = "tab:red"

zorders = np.vectorize(lambda dataset: tp.dataset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% theory
phi_m = 0.585
phi_plot = np.logspace(-3, phi_m, 500)
eta = Krieger_viscosity(phi_plot, phi_m)
Fr_th0 = Froude(0, eta, 7e4, Re_c=350)

# %% masks for plot
# alphas = [0, 45]
alpha_pad = 1.5
alpha0 = 0
mask = ((alpha > alpha0 - alpha_pad) & (alpha < alpha0 + alpha_pad))[plot_idxs]

for ifig in range(3):
    figsize = (1.3 * quarto.regular_fig_width, quarto.regular_fig_height)
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize, sharex=True)
    # inset ax
    axins = ax.inset_axes([0.43, 0.53, 0.55, 0.45])
    axins.set_ylim(0, 0.8)
    axins.set_xlim(right=0.8)
    axins.set_xlabel(r"$\phi$", labelpad=0)
    axins.set_ylabel(r"$\mathcal{F}r$", labelpad=2)
    tp.mscatter(
        phi[plot_idxs][mask],
        Fr[plot_idxs][mask],
        ax=axins,
        m=markers[plot_idxs][mask],
        facecolors=facecolors[plot_idxs][mask],
        edgecolors=edgecolors[plot_idxs][mask],
        lw=0.5,
    )
    #
    tp.mscatter(
        phi[plot_idxs][mask],
        Fr[plot_idxs][mask],
        ax=ax,
        m=markers[plot_idxs][mask],
        facecolors=facecolors[plot_idxs][mask],
        edgecolors=edgecolors[plot_idxs][mask],
        lw=0.5,
    )
    #
    ax.set_ylabel(r"Slumping velocity, $\mathcal{F}r$")
    ax.set_ylim([0, 1.59])

    if ifig > 1:
        ax.plot(phi_plot, Fr_th0, ls="-", color="w")
        axins.plot(phi_plot, Fr_th0, ls="-", color="w")

    if ifig > 0:
        ax.axvline(tp.phi_c, ls=":", lw=1, color="w")
        axins.axvline(tp.phi_c, ls=":", lw=1, color="w")

    ax.set_xlabel(r"Volume fraction, $\phi$")
    ax.set_xscale("log")
    ax.set_xlim([0.0035, 1.2])
    #
    leg1 = fig.legend(
        handles=tp.legend_datasets,
        ncol=1,
        title="Datasets",
        loc="outside right upper",
    )
    leg2 = fig.legend(handles=tp.legend_particles, ncol=1, title="Particles", loc="outside right center")

    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
