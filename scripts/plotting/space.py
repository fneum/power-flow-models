"""
Plotting functions to portray feasible space in loss/flow domain.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D


def plot_feasible_space(line, fn=None):

    print(line[["s_nom", "s_max_pu", "r_pu_eff", "s_nom_max"]])

    fig, ax = plt.subplots(figsize=(5, 4.5))

    x = np.linspace(-line.s_nom, line.s_nom, 100)

    plt.plot(x, line.r_pu_eff * x ** 2, c="k", label=r"$\psi=rp^2$")

    plt.axhline(0, c="firebrick", label=r"$\psi\geq 0$")

    max_loss = line.r_pu_eff * (line.s_max_pu * line.s_nom) ** 2
    ax.axhline(max_loss, c="darkseagreen", label=r"$\psi\leq r(\bar{p}P)^2$")

    flow_upper = -x + line.s_max_pu * line.s_nom
    flow_lower = x + line.s_max_pu * line.s_nom
    plt.plot(x, flow_upper, c="navy", linestyle="--", label=r"$\psi+|p| \leq \bar{p}P$")
    plt.plot(x, flow_lower, c="navy", linestyle="--")

    tangents = []
    for k in [1, 2]:
        p_k = k / 2 * line.s_max_pu * line.s_nom
        loss_k = line.r_pu_eff * p_k ** 2
        slope_k = 2 * line.r_pu_eff * p_k
        offset_k = loss_k - slope_k * p_k
        for sign in [-1, 1]:
            tangent = sign * slope_k * x + offset_k
            tangents.append(tangent)
            plt.plot(x, tangent, c="k", linestyle=":")
    min_loss = [min(max_loss, max(0, max(i))) for i in zip(*tangents)]
    plt.fill_between(
        x, max_loss, min_loss, alpha=0.2, color="firebrick", label="Feasible Space"
    )

    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color="k", linestyle=":", label="Tangents"))

    plt.xlabel("Line Flow [MW]")
    plt.ylabel("Line Losses [MW]")
    plt.ylim([-5, max_loss / line.s_max_pu ** 2 * 0.8])
    plt.legend(handles=handles)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
