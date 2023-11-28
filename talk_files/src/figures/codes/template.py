import cmocean as cmo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# %% colors

# cmap_images = cmo.cm.ice
cmap_images = cmo.cm.gray
cmap_images2 = cmo.cm.gray_r
# cmap_images = cmo.cm.dense_r
# cmap_images = cmr.arctic

color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# color_setups = {'Julien': '#823329', 'Jean': '#FE7F2D',
#                 'Cyril': '#FCCA46', 'Cyril/Marie': '#619B8A', 'Rastello': '#A1C181'}

datasets = {
    "Julien Chauchat": "SedFoam",
    "Jean Schneider": "LEMTA",
    "Cyril Gadal": "IMFT",
    "Marie Rastello - Cyril Gadal": "LEGI",
    "Marie Rastello": "LEGI",
}

# color_datasets = {'SedFoam': '#08415C', '3': '#FE7F2D',
#                   '1': '#FCCA46', '2': '#83B692'}

# color_datasets = {"4": "#780116", "3": "#FE7F2D", "1": "#FCCA46", "2": "#83B692"}
color_datasets = {"SedFoam": "#780116", "LEMTA": "#FE7F2D", "IMFT": "#FCCA46", "LEGI": "#83B692"}


marker_style = {
    "glass beads": "o",
    "silica sand": "h",
    "Hydrogels": "*",
    "PMMA": "D",
    "polystyren beads": "X",
    "SedFoam": "s",
}

legend_names = {
    "glass beads": "glass beads",
    "silica sand": "silica sand",
    "Hydrogels": "hydrogels",
    "PMMA": "PMMA",
    "polystyren beads": "polystyren beads",
    "SedFoam": "simulations",
}

dataset_zorder = {"IMFT": 0, "LEGI": 1, "LEMTA": 2, "SedFoam": 2}
dataset_zorder2 = {"IMFT": 0, "LEGI": 2, "LEMTA": 2, "SedFoam": 1}

# %% colors for other study
cmap_slope = cmo.cm.ice
cmap_slope2 = cmo.cm.algae
color_list = [
    "#006BA4",
    "#FF800E",
    "#ABABAB",
    "#595959",
    "#5F9ED1",
    "#C85200",
    "#898989",
    "#A2C8EC",
    "#FFBC79",
    "#CFCFCF",
]
colors = {
    # #### slope
    "Sand120m_Theta0": cmap_slope(0.9),
    "Slope1": cmap_slope(0.8),
    "Slope3": cmap_slope(0.65),
    "Slope5": cmap_slope(0.45),
    "sand80m_H19": cmap_slope(0.3),
    # 'sand80m_H19': 'tab:cyan',
    "Theta7": cmap_slope2(0.8),
    "Theta10": cmap_slope2(0.5),
    "Theta15": cmap_slope2(0.3),
    # #### settling velocity
    "Saline": color_list[5],
    "Silibeads40_70": color_list[1],
    "silibeads40_70m": color_list[1],
    "Silibeads100_200": color_list[2],
    "Silibeads150_250": color_list[-2],
    # 'silibeads200m_300m': color_list[7]
    "silibeads200m_300m": "tab:purple",
}

# %% corresponding legend

legend_datasets = [
    Line2D(
        [0], [0], marker="s" if dataset == "SedFoam" else "o", color=color_datasets[dataset], ls="none", label=dataset
    )
    for dataset in sorted(color_datasets.keys())
]

legend_particles = [
    Line2D(
        [0],
        [0],
        marker=marker_style[particle],
        markerfacecolor="none",
        markeredgecolor="w",
        ls="none",
        label=legend_names[particle],
    )
    for particle in sorted(marker_style.keys())
]

# %% variables

phi_c = 0.45


# %% plot functions


def mscatter(x, y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers

    if not ax:
        ax = plt.gca()
    sc = ax.scatter(x, y, **kw)
    if (m is not None) and (len(m) == len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc
