"""
Functions yielding statistics about solved networks.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import pandas as pd

from .utils import aggregate_costs


def check_curtailment(n):
    possible = (
        (n.generators_t.p_max_pu * n.generators.p_nom_opt)
        .filter(regex="(ror|solar|wind)", axis=1)
        .sum()
        .sum()
    )
    generated = n.generators_t.p.filter(regex="(ror|solar|wind)", axis=1).sum().sum()
    return (1 - generated / possible) * 100  # %


def check_energy_balance(n):
    gen = n.generators_t.p.multiply(n.snapshot_weightings, axis=0).sum().sum()
    sto = n.storage_units_t.p.multiply(n.snapshot_weightings, axis=0).sum().sum()
    load = n.loads_t.p_set.multiply(n.snapshot_weightings, axis=0).sum().sum()
    if "loss" in n.lines_t.keys():
        line_loss = n.lines_t.loss.multiply(n.snapshot_weightings, axis=0).sum().sum()
    else:
        line_loss = 0
    link_loss = (
        (n.links_t.p0 + n.links_t.p1)
        .multiply(n.snapshot_weightings, axis=0)
        .sum()
        .sum()
    )
    loss = link_loss + line_loss
    loss_share = loss / (gen + sto) * 100  # %
    balance = (gen + sto - loss - load) / 1e6  # TWh
    return balance, loss_share, link_loss / 1e6, line_loss / 1e6


def check_energy_transmitted(n, branch="lines"):
    branches_t = getattr(n, branch + "_t")
    branches = getattr(n, branch)
    return (
        branches_t.p0.abs()
        .multiply(n.snapshot_weightings, axis=0)
        .multiply(branches.length)
        .sum()
        .sum()
        / 1e12
    )  # EWhkm


def check_energy_generated(n):
    n.generators["energy"] = n.generators_t.p.multiply(
        n.snapshot_weightings, axis=0
    ).sum()
    generated = n.generators.groupby("carrier").energy.sum() / 1e6  # TWh
    generated["inflow"] = (
        n.storage_units_t.inflow.multiply(n.snapshot_weightings, axis=0).sum().sum()
        / 1e6
    )  # TWh
    return generated


def check_costs(n):

    abs_c = aggregate_costs(n).sum()  # bn EUR/a
    rel_c = (
        abs_c
        / n.loads_t.p_set.multiply(n.snapshot_weightings, axis=0).sum().sum()
        * 1e9
    )  # EUR/MWh

    return abs_c, rel_c


def check_capacities(n):

    capacities = pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            n.storage_units.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            pd.Series(
                {"links": n.links.eval("length * (p_nom_opt - p_nom) / 1e6 / 2").sum()}
            ),  # TWkm  # links are split and p_nom_opt thus counted double
            pd.Series(
                {"lines": n.lines.eval("length * (s_nom_opt - s_nom) / 1e6").sum()}
            ),  # TWkm
        ]
    )

    if "load" in capacities.index:
        capacities.drop("load", inplace=True)

    return capacities


def check_flow_errors(n, n_pf):

    pf = n_pf.lines_t.p0.stack()

    lopf = n.lines_t.p0.stack()

    mse = pf.sub(lopf).pow(2).mean()
    rmse = mse ** 0.5
    mape = pf.sub(lopf).div(pf).abs().mean()
    mae = pf.sub(lopf).abs().mean()
    corr = pf.corr(lopf)
    r2 = corr ** 2
    return (rmse, mae, mape, corr, r2)


def check_slack(n_pf, logs):
    converged = logs.converged.all(axis=1)

    valley_q = n_pf.generators_t.q.loc[converged].sum(axis=1).min() / 1e3  # GW
    peak_q = n_pf.generators_t.q.loc[converged].sum(axis=1).max() / 1e3  # GW
    tvarh = (
        n_pf.generators_t.q.loc[converged]
        .multiply(n_pf.snapshot_weightings, axis=0)
        .sum()
        .sum()
        / 1e6
    )  # Tvarh

    slack = (n_pf.generators_t.p - n_pf.generators_t.p_set).loc[converged]
    twh = slack.multiply(n_pf.snapshot_weightings, axis=0).sum().sum() / 1e6  # TWh
    peak_p = slack.sum(axis=1).max() / 1e3  # GW
    valley_p = slack.sum(axis=1).min() / 1e3  # GW

    return (twh, peak_p, valley_p, tvarh, peak_q, valley_q)
