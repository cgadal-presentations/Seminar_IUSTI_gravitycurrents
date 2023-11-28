import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg


def x_ad(t_ad, Fr, tau, epsilon):
    return Fr * (t_ad - (1 / tau) * t_ad**2 + epsilon * t_ad**3)


t_ad = np.linspace(0, 23, 300)
xad = x_ad(t_ad, 0.5, 120, epsilon=-4e-4)

figsize = (0.55 * quarto.regular_fig_width, 0.75 * quarto.regular_fig_height)


fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)

ax.plot(t_ad, xad)

ax.axvline(9)
ax.hlines(0.6, xmin=1.27, xmax=3)
ax.vlines(3, ymin=0.6, ymax=1.45)
ax.text(3.5, 0.9, r"$\mathcal{F}_{r} = u_{\rm c}/u_{0}$")

ax.annotate(
    "?",
    xy=(7, 6),
    xytext=(11, 6),
    arrowprops=dict(facecolor="w", arrowstyle="<-", lw=2, shrinkA=0, shrinkB=0),
    ha="center",
    va="center",
)
# ax.text(11, 6, "?", )

ax.set_ylabel(r"Front position, $x_{\rm f}/l_{0}$")
ax.set_xlabel(r"Time, $t/t_{0}$")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.set_xticks([9])
ax.set_xticklabels([r"$\sim \tau$"])
ax.set_yticks([])
ax.spines[["top", "right"]].set_visible(False)

# %% saving figure
figname = "../{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""))
fig.savefig(figname, dpi=600)

fig = sg.fromfile(figname)
newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

fig.set_size(newsize)
fig.save(figname)
