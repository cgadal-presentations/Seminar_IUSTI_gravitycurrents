import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from netCDF4 import Dataset

# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs])

alpha, phi = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
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

mask = (authors == "Marie Rastello") & (particles == "PMMA") & (phi < 0.02)

slopes = [0, 24.0, 31.9, 45.4]
inds = [-1, 1, 1, 1]


figsize = (0.6 * quarto.regular_fig_width, quarto.regular_fig_height)
fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

for slope, ind in zip(slopes, inds):
    id = np.argwhere(mask & (alpha == slope)).squeeze()[ind]
    #
    t = datasets[id]["t"][:].data
    x = datasets[id]["x_front"][:].data
    l0 = datasets[id]["L0"][:].data
    t0 = datasets[id]["t0"][:].data
    #
    print(datasets[id]["phi"][:].data)

    ax.plot(t / t0, x / l0, label=f"{slope:.0f}")

ax.set(xlabel=r"Time, $t/t_{0}$", ylabel=r"Front position, $x_{\rm f}/l_{0}$")
ax.legend(title=r"$\alpha [^\circ]$")

fig.align_labels()

figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=600)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
