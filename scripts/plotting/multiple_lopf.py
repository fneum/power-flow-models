"""
Plotting functions regarding multiple LOPF networks at once.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from .utils import aggregate_costs, assign_carriers


def process_logs(logs):
    attrs = ["time", "peak_mem"]
    df = pd.DataFrame(columns=logs.keys(), index=attrs)
    for fm, mem in logs.items():
        df.at["time", fm] = mem.index[-1] / 3600  # h
        df.at["peak_mem", fm] = mem.max() / 1e3  # GB
    return df


def plot_performance(logs, attr, model_names=None, colors="forestgreen", fn=None):

    df = process_logs(logs)

    if model_names is not None:
        df.rename(columns=model_names, inplace=True)

    fig, ax = plt.subplots(figsize=(4.5, 2.5))

    df.T[attr].plot.bar(ax=ax, color=colors)

    if attr == "peak_mem":
        plt.ylabel("Peak Memory [GB]")
    else:
        plt.ylabel("Solving Time [h]")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_cost_bar(networks, model_names, fn=None):

    for n in networks.values():
        assign_carriers(n)

    costs = pd.concat({k: aggregate_costs(v) for k, v in networks.items()}, axis=1).T

    if "load" in costs.columns:
        costs.drop(columns=["load"], inplace=True)

    colors = n.carriers.color.reindex(index=costs.columns).values

    costs.rename(columns=n.carriers.nice_name, index=model_names, inplace=True)
    costs.columns.name = "Technology"

    fig, ax = plt.subplots(figsize=(8, 4))

    costs.plot.bar(ax=ax, stacked=True, color=colors)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], ncol=1, bbox_to_anchor=(1, 1.2))

    plt.xticks(rotation=0)
    plt.ylabel("Total System Costs [bn Euro / a]")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def optimised_capacities(n, c, regex="()"):
    attr = "s" if c == "Line" else "p"
    return n.df(c)[f"{attr}_nom_opt"].filter(regex=regex)


def plot_capacity_correlation(
    networks, c, model_names, regex="", triangle=False, fn=None
):

    regexb = "(" + regex + ")"

    df = pd.DataFrame(
        {fm: optimised_capacities(n, c, regexb) for fm, n in networks.items()}
    )
    df.rename(columns=model_names, inplace=True)

    corr = df.corr()

    mask = None
    if triangle:
        mask = np.triu(np.ones_like(corr, dtype=np.bool))

    fig, ax = plt.subplots(figsize=(4, 4))

    sns.heatmap(
        df.corr(),
        vmin=0.5,
        mask=mask,
        cmap="viridis",
        square=True,
        annot=True,
        fmt=".2",
        ax=ax,
        cbar=False,
    )

    auxn = next(iter(networks.values()))

    if regex != "":
        plt.title(getattr(auxn.carriers.nice_name, regex, regex))
    else:
        plt.title(c)

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_price_duration_curve(networks, model_names, ignore=[], fn=None):

    fig, axs = plt.subplots(
        1, 2, sharey=True, figsize=(10, 3), gridspec_kw={"width_ratios": [1, 5]}
    )

    for k, v in networks.items():
        if k in ignore:
            continue
        y = v.buses_t.marginal_price.stack().sort_values().reset_index(drop=True)
        y.index = [100 * i / len(y) for i in y.index]
        for ax in axs:
            y.plot(label=model_names[k], ax=ax)
    
    axs[0].set_xlim([-0.01, 0.2])
    axs[1].set_xlim([80, 100])
    
    axs[0].set_ylabel("Nodal price [EUR/MWh]")
    axs[1].set_xlabel("Share of Snapshots and Nodes [%]")

    plt.tight_layout()
    plt.legend(loc="upper left", ncol=2)
    plt.ylim([-100, 1000])

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
