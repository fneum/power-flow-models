"""
Plotting functions regarding a single network using both LOPF and PF outputs.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .utils import reference, plot_hist_helper, line_loading


def plot_loss_comparison(n, n_pf, style="hist2d", title="", fn=None, norm="max"):

    fig, ax = plt.subplots(figsize=(6, 5))

    max_loss = n.lines.r_pu_eff * (n.lines.s_max_pu * n.lines[f"s_nom_{norm}"]) ** 2

    pf_loss = ((n_pf.lines_t.p0 + n_pf.lines_t.p1) / max_loss).stack()
    lopf_loss = (n.lines_t.loss / max_loss).stack()

    xlim = [0, 1]
    ylim = [0, 1.2]

    plot_hist_helper(ax, lopf_loss, pf_loss, xlim, ylim, style=style, vmax=100)

    reference(ax, *xlim, f=lambda x: x)

    if style in ["hexbin", "hist2d"]:
        cb = plt.colorbar(ax=ax, shrink=0.95)
        cb.set_label("Count")

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    plt.xlabel("Rel. Losses (LOPF)")
    plt.ylabel("Rel. Losses (PF)")

    plt.title(title)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_flow_comparison(n, n_pf, style="hist2d", title="", fn=None):

    pf = n_pf.lines_t.relative_loading.stack()
    lopf = n.lines_t.relative_loading.stack()

    fig, ax = plt.subplots(figsize=(6, 5))

    xlim = [-1.1, 1.1]
    ylim = [-0.8, 0.8]

    plot_hist_helper(ax, pf, lopf, xlim, ylim, style=style, vmax=600)

    if style in ["hexbin", "hist2d"]:
        cb = plt.colorbar(ax=ax, shrink=0.75)
        cb.set_label("Count")

    reference(ax, -1, 1)

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    ax.set_yticks(np.arange(ylim[0], ylim[1] + 0.1, 0.2))
    ax.set_xticks(np.arange(xlim[0] + 0.1, xlim[1] + 0.1, 0.2))
    plt.xticks(rotation=90)

    ax.set_aspect(1)

    plt.xlabel("Rel. Line Flows (PF)")
    plt.ylabel("Rel. Line Flows (LOPF)")

    plt.title(title)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_duration_curve(n, n_pf, apparent=True, fn=None):

    fig, ax = plt.subplots(figsize=(5, 4))

    def duration_curve(nc, s=True, label=""):
        series = pd.Series(
            line_loading(nc, apparent=s).stack().sort_values(ascending=False).values
        )
        series.index = [i / len(series) * 100 for i in series.index]
        series.plot(ax=ax, label=label)

    duration_curve(n, s=False, label="LOPF (Approx. power flow)")
    duration_curve(n_pf, s=apparent, label="PF (AC power flow)")

    plt.ylim([-0.1, 1.4])
    plt.legend(loc="upper right")
    plt.ylabel("Relative Line Loading [-]")
    plt.xlabel("Share of Snapshots [%]")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
