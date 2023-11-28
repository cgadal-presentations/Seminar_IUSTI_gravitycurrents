import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
import template as tp
from lmfit import Model
from matplotlib.colors import to_rgba
from netCDF4 import Dataset

plt.rcParams["figure.constrained_layout.hspace"] = 0
plt.rcParams["figure.constrained_layout.h_pad"] = 0.0005


def lambda_var(x, a, c, th):
    return np.piecewise(x, [x < th, x >= th], [lambda x: c, lambda x: a * (x - th) + c])


# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, St, H0, L0, Fr, lamb = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
            d.variables["St"][:].data,
            d.variables["H0"][:].data,
            d.variables["L0"][:].data,
            d.variables["Fr"][:].data,
            d.variables["lamb"][:].data,
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
a = H0 / L0

# %% graphic specifications
dataset_idx = np.vectorize(tp.datasets.get)(authors)

mask_phi = phi < tp.phi_c
alphas = np.ones_like(Fr)
alphas[~mask_phi] = 0.4

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == "SedFoam"] = "s"

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
facecolors = np.array([to_rgba(c, a) for c, a in zip(facecolors, alphas)])

mask_nosuspended = (authors == "Marie Rastello") & (H0 / Ha < 1)
edgecolors = np.array([to_rgba("dimgray", a) for a in alphas])
edgecolors[mask_nosuspended] = np.array([to_rgba("tab:red", 0.4) for a in alphas[mask_nosuspended]])

zorders = np.vectorize(lambda dataset: tp.dataset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# alpha0 = [0, 7, 15]
alpha_pad = 1.5
a0 = 7

# %% fit for alpha = 7

model = Model(lambda_var)
params = model.make_params()
p0 = {"a": 1, "c": 0, "th": 1e-2}

for par in params.keys():
    params[par].set(value=p0[par])

params["th"].vary = False
params["c"].vary = False

mask_alpha = (alpha > a0 - alpha_pad) & (alpha < a0 + alpha_pad)

result = model.fit(lamb[mask_alpha], params, x=(St)[mask_alpha])

# # %% Figure
for ifig in range(3):
    figsize = (1.3 * quarto.regular_fig_width, quarto.regular_fig_height)
    fig, axarr = plt.subplots(2, 1, constrained_layout=True, figsize=figsize)
    mask_alpha = (alpha > a0 - alpha_pad) & (alpha < a0 + alpha_pad)
    mask = (mask_alpha)[plot_idxs]
    #
    tp.mscatter(
        (St)[plot_idxs][mask],
        Fr[plot_idxs][mask],
        ax=axarr[0],
        m=markers[plot_idxs][mask],
        facecolors=facecolors[plot_idxs][mask],
        edgecolors=edgecolors[plot_idxs][mask],
        lw=0.5,
    )
    #
    moy, std = np.nanmean(Fr[mask_alpha & mask_phi]), np.nanstd(Fr[mask_alpha & mask_phi])
    axarr[0].axhline(moy, color="w", ls=":", zorder=-10)
    if ifig > 0:
        tp.mscatter(
            (St)[plot_idxs][mask],
            lamb[plot_idxs][mask],
            ax=axarr[1],
            m=markers[plot_idxs][mask],
            facecolors=facecolors[plot_idxs][mask],
            edgecolors=edgecolors[plot_idxs][mask],
            lw=0.5,
        )
        axarr[1].axhline(0, color="w", ls=":", zorder=-10)
        if ifig > 1:
            x_plot = np.logspace(np.log10(result.params["th"]), 0, 300)
            axarr[1].plot(
                x_plot,
                result.params["a"] * (x_plot - result.params["th"]) + result.params["c"],
                color="w",
                ls="-",
                label=r"$1/\tau \propto \mathcal{S}t$",
            )
            axarr[1].legend(loc="upper left")
    else:
        axarr[1].set_axis_off()

    axarr[0].set_ylabel(r"Slumping velocity, $\mathcal{F}r$")
    axarr[1].set_ylabel(r"Attenuation, $1/\tau$")

    axarr[1].set_xscale("log")
    axarr[0].set_xscale("log")

    axarr[0].set_ylim(0, 1.59)
    axarr[1].set_ylim(-0.02, 0.07)
    axarr[1].set_xlim(3e-4, 1)
    axarr[0].set_xlim(3e-4, 1)
    #
    axarr[1].set_yticks([0, 0.03, 0.06])
    #
    axarr[1].set_xlabel(r"Stokes number, $\mathcal{S}t$")
    if ifig == 0:
        axarr[0].set_xlabel(r"Stokes number, $\mathcal{S}t$")
    else:
        axarr[0].set_xticks([])
    #
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
