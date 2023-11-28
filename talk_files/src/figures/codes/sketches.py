import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import PyThemes.quarto as quarto
import svgutils.transform as sg
from matplotlib.path import Path


def sind(x):
    return np.sin(x * np.pi / 180.0)


def cosd(x):
    return np.cos(x * np.pi / 180.0)


def Rotation_matrix(theta):
    return np.array([[cosd(theta), -sind(theta)], [sind(theta), cosd(theta)]])


# %%
# ## Sketches parameters
# colors
# color_water = 'aliceblue'
color_water = "lightcyan"
# color_water = "#8df2f2"
# color_water_salt = 'lightgreen'
color_water_salt = "#F3F7D4"
color_sed = "peru"
color_walls = "k"
alpha_water = 1
color_mixing = "grey"
color_text = "k"

# dimension parameters
tank_height = 45  # cm
tank_length = 165  # cm
door_pos = 15  # cm
door_height = 1.1 * tank_height
door_pad = 0 * tank_height
y_bottom = 0  # cm, at the end of the canal
x_bottom = 0  # cm, at the end of the canal
mixing_height = tank_height
mixing_pad = 0.2 * tank_height
mixing_width = 0.15 * mixing_height


SLOPES = np.array([-7, -40, 0, -7])
water_height_facts = [0.89, 1, 0.89, 0.89]  # proportion of tank height
width_ratios = tank_length * cosd(SLOPES) + np.abs(tank_height * sind(SLOPES))

figsizes = [
    (0.6 * quarto.regular_fig_width, 0.5 * quarto.regular_fig_height),
    (0.6 * quarto.regular_fig_width, 0.75 * quarto.regular_fig_height),
    (0.6 * quarto.regular_fig_width, 0.5 * quarto.regular_fig_height),
    (0.75 * quarto.regular_fig_width, 0.5 * quarto.regular_fig_height),
]
labels = ["IMFT", "LEGI", "LEMTA", "general"]
path_logos = "../../images/logos"
logo_paths = [
    os.path.join(path_logos, "logo_IMFT.png"),
    os.path.join(path_logos, "logo_LEGI.png"),
    os.path.join(path_logos, "logo_LEMTA.png"),
]
logos = [plt.imread(logo) for logo in logo_paths] + [None]
bounds = [[0.65, 0.7, 0.35, 0.3], [0.65, 0.8, 0.35, 0.2], [0.65, 0.75, 0.35, 0.25], None]

# %%%%%%% Sketches

for ifig, (figsize, label, slope, water_height_fact, logo, b) in enumerate(
    zip(figsizes, labels, SLOPES, water_height_facts, logos, bounds)
):
    fig, ax = plt.subplots(1, 1, layout="constrained", figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    #
    ax.set_aspect("equal")
    ax.set_axis_off()
    # ax.set_xticks([])
    # ax.set_yticks([])
    if b is not None:
        ax_logo = ax.inset_axes(b)
        ax_logo.set_axis_off()
        ax_logo.imshow(logo)

    # parameters
    water_height = water_height_fact * tank_height  # cm, at end of the canal

    # ## reference points
    slope_vec = np.array([cosd(slope), sind(slope)])
    slope_vec_up = np.array([-sind(slope), cosd(slope)])
    down_vec = np.array([0, -1])
    # tank
    bottom_right = np.array([x_bottom, y_bottom])
    bottom_left = bottom_right - tank_length * slope_vec
    top_right = bottom_right + tank_height * slope_vec_up
    top_left = bottom_left + tank_height * slope_vec_up
    # door
    bottom_door = bottom_left + door_pos * slope_vec + door_pad * slope_vec_up
    top_door = bottom_door + door_height * slope_vec_up
    # water (ambient)
    bottom_door_water = bottom_door
    bottom_right_water = bottom_right
    top_right_water = bottom_right_water + water_height * slope_vec_up
    if ifig == 1:
        top_door_water = top_right_water - (tank_length - door_pos) * slope_vec
    else:
        top_door_water = top_right_water - (tank_length - door_pos) * np.array([1, 0]) / cosd(slope)
    #
    # water (reservoir)
    bottom_right_water_reservoir = bottom_door_water
    top_right_water_reservoir = top_door_water
    bottom_left_water_resevoir = bottom_left
    top_left_water_reservoir = top_door_water - door_pos * np.array([1, 0]) / cosd(slope)

    xy_water = [bottom_door_water, top_door_water, top_right_water, bottom_right_water]
    xy_water_reservoir = [
        bottom_left_water_resevoir,
        top_left_water_reservoir,
        top_right_water_reservoir,
        bottom_right_water_reservoir,
    ]

    # mixing
    bottom_mixing = bottom_left + 0.5 * door_pos * slope_vec + mixing_pad * slope_vec_up
    top_mixing = bottom_mixing + mixing_height * slope_vec_up
    bottom_left_mixing = bottom_mixing - 0.5 * mixing_width * slope_vec
    bottom_right_mixing = bottom_mixing + 0.5 * mixing_width * slope_vec

    # sediment position generation
    ngrains = 70
    # np.random.seed(220212021)
    np.random.seed(999999)
    xsed = door_pos * np.random.random((ngrains,))
    ysed = 1.2 * water_height * np.random.random((ngrains,))
    # if ifig == 0:
    #     xsed, ysed = door_pos * \
    #         np.random.random((ngrains, )), water_height * \
    #         (1+sind(slope))*np.random.random((ngrains, ))
    # else:
    #     xsed, ysed = door_pos*np.random.random((ngrains, )), 0.5*water_height*(
    #         1+sind(slope))*np.random.random((ngrains, ))
    xsed, ysed = (np.dot(Rotation_matrix(slope), np.array([xsed, ysed])).T - tank_length * slope_vec).T

    # ## ploting sketch
    hpad = 0.007 * tank_length
    ax.set_xlim(bottom_left[0] - 0.2 * door_height * np.abs(sind(slope)) - hpad, top_right[0] + hpad)
    ax.set_ylim(bottom_right[1] - 0.2 * door_height, top_door[1] + 0.4 * door_height)
    # ax.set_xlim(bottom_left[0] - 0.05*tank_length,
    #             top_right[0] + 0.05*tank_length)
    # ax.set_ylim(bottom_right[1] - 0.2*door_height, top_door[1]+0.4*door_height)

    #
    # ## tank walls
    ax.plot([bottom_left[0], bottom_right[0]], [bottom_left[1], bottom_right[1]], color=color_walls)
    ax.plot([bottom_left[0], top_left[0]], [bottom_left[1], top_left[1]], color=color_walls)
    ax.plot([bottom_right[0], top_right[0]], [bottom_right[1], top_right[1]], color=color_walls)
    ax.plot([bottom_door[0], top_door[0]], [bottom_door[1], top_door[1]], color=color_walls)
    if ifig == 1:
        ax.plot([top_door_water[0], top_right[0]], [top_door_water[1], top_right[1]], color=color_walls)
    # ## water
    poly_water = plt.Polygon(xy_water, facecolor=color_water, alpha=alpha_water, edgecolor=None)
    ax.add_patch(poly_water)
    poly_water_reservoir = plt.Polygon(
        xy_water_reservoir,
        facecolor=color_water_salt if ifig >= 2 else color_water,
        alpha=alpha_water,
        edgecolor=None,
    )
    ax.add_patch(poly_water_reservoir)

    # ## sediments
    mask_in = Path(xy_water_reservoir).contains_points(np.array([xsed, ysed]).T, radius=3)
    ax.scatter(xsed[mask_in], ysed[mask_in], color=color_sed, s=0.7)
    # ax.scatter(xsed[ysed < water_height], ysed[ysed <
    #            water_height], color=color_sed, s=0.7)

    # # ## mixing
    # ax.plot([bottom_mixing[0], top_mixing[0]], [
    #         bottom_mixing[1], top_mixing[1]], color=color_mixing)
    # ax.plot([bottom_left_mixing[0], bottom_right_mixing[0]],
    #         [bottom_left_mixing[1], bottom_right_mixing[1]], color=color_mixing, lw=2)

    # ## annotations
    ax.plot(
        [bottom_right[0], bottom_right[0] - 0.4 * tank_length], [bottom_right[1], bottom_right[1]], ls="--", color="k"
    )
    theta = np.linspace(180, 180 + slope, 100)
    x, y = 0.35 * tank_length * np.array([cosd(theta), sind(theta)]) + bottom_right[:, None]
    ax.plot(x, y, color="k")
    #
    if slope != 0:
        ax.text(
            bottom_right[0] - 0.4 * tank_length, 0.75 * y.mean(), r"$\alpha$", ha="right", va="center", color=color_text
        )

    ax.annotate(
        "",
        xytext=top_door,
        xy=top_door + 0.4 * door_height * slope_vec_up,
        arrowprops=dict(arrowstyle="-|>", shrinkA=4, shrinkB=0, color="k"),
    )

    # xy = bottom_right - 0.195 * door_height * slope_vec_up
    # xytext = bottom_door - 0.195 * door_height * slope_vec_up
    # ax.annotate("", xytext=xytext, xy=xy, arrowprops=dict(arrowstyle="<->", shrinkA=0, shrinkB=0, color="k"))
    # xytext = (xy + xytext) / 2 - slope_vec_up * 0.05 * door_height
    # ax.text(xytext[0], xytext[1], r"$L_{1}$", ha="center", va="top", color=color_text)
    #
    xy = bottom_left - 0.195 * door_height * slope_vec_up
    xytext = bottom_door - 0.195 * door_height * slope_vec_up
    ax.annotate("", xytext=xytext, xy=xy, arrowprops=dict(arrowstyle="<->", shrinkA=0, shrinkB=0, color="k"))
    xytext = (xy + xytext) / 2
    ax.text(xytext[0], xytext[1] - 0.06 * door_height, r"$l_{0}$", ha="center", va="top", color=color_text)
    #
    xy = bottom_door + 0.1 * door_height * slope_vec
    xytext = top_right_water_reservoir + 0.1 * door_height * slope_vec
    ax.annotate("", xytext=xytext, xy=xy, arrowprops=dict(arrowstyle="<->", shrinkA=0, shrinkB=0, color="k"))
    xytext = (xy + xytext) / 2 + 0.1 * door_height * slope_vec
    ax.text(xytext[0], xytext[1], r"$h_{0}$", ha="center", va="center", color=color_text)

    figname = "../{}_{}.svg".format(sys.argv[0].split(os.sep)[-1].replace(".py", ""), label)
    fig.savefig(figname, dpi=600, facecolor="white")

    fig = sg.fromfile(figname)
    newsize = ["{}pt".format(quarto.scaling_factor * float(i.replace("pt", ""))) for i in fig.get_size()]

    fig.set_size(newsize)
    fig.save(figname)
