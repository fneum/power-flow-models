"""
Utility functions for plotting functionality.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_memory(fn):
    memlog = pd.read_csv(
        fn, sep=" ", index_col=1, usecols=[1, 2], header=None
    ).squeeze()
    memlog.index = [t - memlog.index[0] for t in memlog.index]
    return memlog


def aggregate_costs(n, existing_only=False, by_carrier=True):

    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
    )

    costs = {}
    for c in n.iterate_components(components.keys()):
        p_nom, p_attr = components[c.name]
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].multiply(n.snapshot_weightings, axis=0).sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs) / 1e9  # bn EUR/a

    if by_carrier:
        costs = costs.groupby(level=2).sum()

    return costs


def assign_carriers(n):

    if "Load" in n.carriers.index:
        n.carriers = n.carriers.drop("Load")

    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"

    if n.links.empty:
        n.links["carrier"] = pd.Series(dtype=str)

    config = {
        "AC": {"color": "rosybrown", "nice_name": "HVAC Line"},
        "DC": {"color": "darkseagreen", "nice_name": "HVDC Link"},
    }

    for c in ["AC", "DC"]:
        if c in n.carriers.index:
            continue
        n.carriers = n.carriers.append(pd.Series(config[c], name=c))


def line_loading(n, apparent=True, relative=True):
    p = np.maximum(n.lines_t.p0.abs(), n.lines_t.p1.abs())

    if apparent:
        q = np.maximum(n.lines_t.q0.abs(), n.lines_t.q1.abs())
    else:
        q = 0

    s = np.sqrt(p ** 2 + q ** 2)

    if relative:
        return s / n.lines.s_nom_opt
    else:
        return s


def plot_hist_helper(ax, x, y, xlim, ylim, style="hexbin", vmin=0, vmax=400):

    if style == "scatter":
        ax.scatter(x, y, alpha=0.2, marker=".")
    elif style == "hexbin":
        scope = xlim + ylim
        plt.hexbin(
            x,
            y,
            gridsize=100,
            cmap=plt.cm.viridis,
            vmin=vmin,
            vmax=vmax,
            extent=xlim + ylim,
        )
    elif style == "hist2d":
        scope = np.array([xlim, ylim])
        plt.hist2d(
            x,
            y,
            bins=100,
            cmap=plt.cm.viridis,
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            range=scope,
        )


def reference(ax, low, high, res=100, lw=1, c="k", f=lambda x: x):
    x = np.linspace(low, high, res)
    ax.plot(x, f(x), c="k", linewidth=lw)
