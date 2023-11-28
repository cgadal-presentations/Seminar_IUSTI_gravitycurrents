import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
import template as tp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset

# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs if Dataset(run).particle_type != "saline water"])

# %% create data vectors
alpha, phi, Re, St, H0, a = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
            d.variables["Re"][:].data,
            d.variables["St"][:].data,
            d.variables["H0"][:].data,
            d.variables["a"][:].data,
        ]
        for d in datasets
    ]
).T

authors, particles = np.array(
    [
        [
            d.author,
            d.particle_type,
        ]
        for d in datasets
    ]
).T

Ha = np.array(
    [d.variables["H_a"][:].data if "H_a" in d.variables.keys() else d.variables["H0"][:].data for d in datasets]
)

# %% graphic vector for plots
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

# %% figure
ylabels = [r"Volume fraction, $\phi$", r"angle, $\alpha~[^\circ]$", r"Stokes number, $\mathcal{S}t$"]

figsize = 1.1 * np.array(quarto.regular_fig_size)
fig, axarr = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)

for ax, var, ylabel in zip(axarr.flatten(), [phi, alpha, St], ylabels):
    tp.mscatter(
        Re[plot_idxs],
        var[plot_idxs],
        ax=ax,
        m=markers[plot_idxs],
        facecolors=facecolors[plot_idxs],
        edgecolors=edgecolors[plot_idxs],
        lw=0.5,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Reynolds number, $\mathcal{R}e$")
    ax.set_ylabel(ylabel)


axarr[0, 1].set_yscale("linear")
axarr[0, 1].set_ylim(bottom=-1)

# legend
divider = make_axes_locatable(axarr[1, 1])
ax_leg2 = divider.append_axes("bottom", size="100%", sharex=axarr[1, 1])
axarr[1, 1].axis("off")
ax_leg2.axis("off")

leg = axarr[1, 1].legend(
    handles=tp.legend_datasets,
    loc="upper center",
    ncol=2,
    borderaxespad=0,
    title="Datasets",
    mode="expand",
)
leg = ax_leg2.legend(
    handles=tp.legend_particles, loc="upper center", ncol=2, borderaxespad=0, title="Particles", mode="expand"
)

fig.align_labels()

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=600)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
