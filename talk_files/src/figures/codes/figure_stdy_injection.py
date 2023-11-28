import os
import sys

import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib import pyplot as plt
from uncertainties import unumpy as unp

path_data = "/media/cyril/IMFT data/202211_VISITE_JEAN/cyril/analysis_images/compute_average_profiles"

# #### parameters
manip_density = [1, 2, 4]
manip_particle = [5, 6, 8]

rho_c = {1: 1078.5, 2: 1024.4, 4: 1049.4, 5: 1047.5, 6: 1076.6, 8: 1017.9, 9: 1081.7}

rho_a = {1: 998.4, 2: 1007.5, 4: 1000.4, 5: 1000.4, 6: 1000.4, 8: 1000.1, 9: 1030.5}

phi = {1: 0, 2: 0, 3: 0, 5: 3, 6: 3, 8: 3, 9: 3}

zref = {
    1: 0.64,
    2: 0.57,
    4: 0.71,
    #
    5: 0.07,
    6: 0.05,
    8: 0.4,
    9: 0.4,
}

rho_p = 1050  # kg/m3
d_p = 500e-6  # m
Q = 1.2  # L/s
W = 21.5  # cm

labels = ["settling", "neutral", "buoyant"]

for ifig in range(2):
    figsize = (1.4 * quarto.regular_fig_width, 0.9 * quarto.regular_fig_height)
    fig, axarr = plt.subplots(1, 3, layout="constrained", sharey=True, figsize=figsize)
    if ifig > 0:
        axarr_phi = np.array([ax.twiny() for ax in axarr])
    for manip_number in manip_density + manip_particle:
        print(manip_number)
        #
        data = np.load(os.path.join(path_data, f"manip{manip_number:0d}", "vertical_profile.npy"), allow_pickle=True)
        drho_1 = rho_c[manip_number] - rho_a[manip_number]
        drho_2 = rho_p - rho_c[manip_number]
        print(drho_1, drho_2)
        #
        if np.abs(drho_1 - 17) < 4:
            axind = 0
        elif np.abs(drho_1 - 50) < 4:
            axind = 1
        elif np.abs(drho_1 - 80) < 4:
            axind = 2
        #
        z, d = unp.nominal_values(data)
        zerr, derr = unp.std_devs(data)
        nsigma = 2
        if manip_number == 4:
            d = d * (rho_a[manip_number] + 50) / (60 + rho_a[manip_number])
        if manip_number in manip_density:
            (a,) = axarr[axind].plot(d - rho_a[manip_number], z - zref[manip_number], color="C1", label="salt water")
            axarr[axind].fill_betweenx(
                z - zref[manip_number],
                d - rho_a[manip_number] - nsigma * derr,
                d - rho_a[manip_number] + nsigma * derr,
                color="C1",
                alpha=0.3,
            )
            ind = np.nanargmin(np.abs(d - rho_p))
            axarr[axind].scatter(d[ind] - rho_a[manip_number], z[ind] - zref[manip_number], color="C1")
        elif ifig > 0:
            (a,) = axarr_phi[axind].plot(d, z - zref[manip_number], color="C2", ls="--", label="particles")
            axarr_phi[axind].fill_betweenx(
                z - zref[manip_number], d - nsigma * derr, d + nsigma * derr, color="C2", alpha=0.3
            )
        #
        axarr[axind].set_xlim(0, drho_1)
        #
    for ax in axarr:
        ax.set_ylim(0, 9)
        ax.set_xlabel(r"$\rho_{\rm c} - \rho_{\rm a}$", color="C1")
        ax.tick_params(axis="x", labelcolor="C1")
        #
    if ifig > 0:
        for ax in axarr_phi:
            ax.set(xlim=(0, 4))
            ax.set_xlabel(r"$\phi~[\%]$", color="C2")
            ax.tick_params(axis="x", labelcolor="C2")
        for ax, title in zip(axarr, labels):
            ax.set_title(title)

    axarr[0].set_ylabel(r"Height, $z$ [cm]")
    axarr[axind].legend()
    #
    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), ifig)
    fig.savefig(figname, dpi=600)

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
