"""

"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

from .utils import line_loading


def plot_network_losses(n, fn=None):

    lc = n.lines_t.loss.mean()
    lw = line_loading(n, apparent=False, relative=False).mean() / 700

    fig, ax = plt.subplots(
        figsize=(7, 7), subplot_kw={"projection": ccrs.PlateCarree()}
    )

    n.plot(
        ax=ax,
        color_geomap=True,
        line_widths=lw,
        line_colors=lc,
        bus_sizes=5e-3,
        bus_colors="darkgray",
        bus_alpha=0.8,
        link_widths=0,
    )

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_v_ang_diff(n_pf, fn=None):

    v_ang_1 = n_pf.buses_t.v_ang.loc[:, n_pf.lines.bus1]
    v_ang_1.columns = n_pf.lines.index

    v_ang_0 = n_pf.buses_t.v_ang.loc[:, n_pf.lines.bus0]
    v_ang_0.columns = n_pf.lines.index

    v_ang_diff = (v_ang_1 - v_ang_0).applymap(lambda x: x * 180 / np.pi)

    fig, ax = plt.subplots(figsize=(5, 3))

    v_ang_diff.stack().plot.hist(bins=np.arange(-90, 90, 5), density=True)

    plt.xticks(np.arange(-90, 90, 10), rotation=90)

    plt.xlabel("Voltage Angle Difference [Degrees]")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")

