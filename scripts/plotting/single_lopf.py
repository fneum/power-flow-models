"""
Plotting functions regarding a single LOPF network.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt
import pandas as pd

from .utils import reference, plot_hist_helper


def plot_flow_vs_loss(n, norm="max", style="hist2d", title="", fn=None):

    fig, ax = plt.subplots(figsize=(6, 5))

    loading = (n.lines_t.p0 / n.lines[f"s_nom_{norm}"] / n.lines.s_max_pu).stack()
    max_loss = n.lines.r_pu_eff * (n.lines.s_max_pu * n.lines[f"s_nom_{norm}"]) ** 2
    relative_loss = (n.lines_t.loss / max_loss).stack()

    xlim = [-1, 1]
    ylim = [0, 1.1]

    plot_hist_helper(ax, loading, relative_loss, xlim, ylim, vmax=100, style=style)

    reference(ax, *xlim, f=lambda x: x ** 2)

    if style in ["hexbin", "hist2d"]:
        cb = plt.colorbar(ax=ax, shrink=0.95)
        cb.set_label("Count")

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    plt.ylabel("Rel. Losses (LOPF)")
    plt.xlabel("Rel. Line Flows (LOPF)")

    plt.title(title)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_negative_marginal_prices(n, fn=None, max_mp=-0.5):

    mp = n.buses_t.marginal_price.stack()

    neg_mp = pd.Series(mp.loc[mp < max_mp].sort_values().values)

    if neg_mp.empty:
        return

    fig, ax = plt.subplots(figsize=(5, 4))
    neg_mp.plot(ax=ax)

    plt.ylabel("EUR/MWh")
    plt.xlabel("Count")
    plt.title(f"Frequency: {len(neg_mp) / len(mp) * 100:f} %")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
