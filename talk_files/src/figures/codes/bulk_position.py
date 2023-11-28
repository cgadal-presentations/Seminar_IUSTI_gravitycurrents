import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
import template as tp
from netCDF4 import Dataset

# %% Load data
path_data = "/home/cyril/Documents/Work/Research/Projects/Turbity_Currents/palagram_monograph/data/output_data"
list_runs = sorted(glob.glob(os.path.join(path_data, "*.nc")))
list_fitresults = sorted(glob.glob(os.path.join(path_data, "fitresult*")))

datasets = np.array([Dataset(run) for run in list_runs])

# %% mask data
particles = np.array([d.particle_type for d in datasets]).T
mask = particles != "saline water"

zorder_setups = {
    "Cyril Gadal": -9,
    "Marie Rastello": -10,
    "Jean Schneider": -6,
    "Julien Chauchat": -7,
    "Marie Rastello - Cyril Gadal": -10,
}

# %% figure

for ifig in range(2):
    figsize = (0.55 * quarto.regular_fig_width, 0.95 * quarto.regular_fig_height)
    fig, axarr = plt.subplots(2, 1, gridspec_kw={"height_ratios": [0.001, 1]}, figsize=figsize, layout="constrained")
    # All runs
    ax = axarr[1]
    if ifig == 0:
        for d in datasets[mask]:
            ax.scatter(
                d.variables["t"][:].data,
                d.variables["x_front"][:].data,
                s=1,
                color=tp.color_datasets[tp.datasets[d.author]],
                marker=tp.marker_style[d.particle_type] if d.author != "Julien" else "s",
                zorder=zorder_setups[d.author],
                rasterized=True,
            )

        ax.set_ylabel(r"Front position, $x_{\rm f}$ [m]")
        ax.set_xlabel(r"Time, $t$ [s]")
        ax.set_xlim(0, 100)
        ax.set_ylim(bottom=0)
    else:
        # selected run, non-dimensional
        i_runs = [0, 158, 94, 197, 210, 202]

        for i in i_runs:
            d = datasets[list_runs.index(os.path.join(path_data, f"run_{i:03d}.nc"))]
            # print(d.author)
            #
            t_ad = d.variables["t0"][:].data
            x_ad = d.variables["L0"][:].data
            #
            x_axis = d.variables["t"][:].data / t_ad
            y_axis = d.variables["x_front"][:].data / x_ad
            ax.plot(x_axis, y_axis, color=tp.color_datasets[tp.datasets[d.author]])
            # ax.scatter(
            #     x_axis,
            #     y_axis,
            #     color=tp.color_datasets[tp.datasets[d.author]],
            #     s=3,
            #     marker=tp.marker_style[d.particle_type] if d.author != "Julien" else "s",
            #     rasterized=True,
            # )
            # # plotting fit
            # fitresult = load_modelresult(os.path.join(path_data, f"fitresult_run_{i:03d}.save"))
            # xplot = np.linspace(fitresult.userkws["t"].min(), fitresult.userkws["t"].max(), 500)
            # ax.plot(xplot, fitresult.eval(t=xplot), color="k", lw=1, ls="--")

        # annotations
        ax.annotate(
            r"",
            xy=(19.3, 10.3),
            xycoords="data",
            xytext=(9, 12),
            textcoords="data",
            arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=0,angleB=120"),
        )
        ax.text(15.5, 12, r"$\alpha \nearrow$")

        ax.annotate(
            r"",
            xy=(29.9, 13.45),
            xycoords="data",
            xytext=(36.6, 8.16),
            textcoords="data",
            arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=-90,angleB=-20"),
        )
        ax.text(31, 10.4, r"$\mathcal{S}t \nearrow$")

        ax.set_ylabel(r"Front position, $x_{\rm f}/l_{0}$")
        ax.set_xlabel(r"Time, $t/t_{0}$")
        ax.set_xlim(-1, 40)
        ax.set_ylim(-0.5, 17)

    axarr[0].axis("off")
    leg = axarr[0].legend(handles=tp.legend_datasets, ncol=2, title="Datasets", loc="upper center", borderaxespad=0)

    fig.align_labels()

    # %% saving figure
    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
