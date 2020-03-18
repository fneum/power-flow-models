"""
Extra functionality implementing loss models.
"""

# TODO implement nomopyomo variant

__author__ = "Fabian Neumann (KIT), Anika Bitsch (KIT)"
__copyright__ = (
    "Copyright 2019-2020 Fabian Neumann (KIT), Anika Bitsch (KIT), GNU GPL 3"
)

import pandas as pd

from pypsa.descriptors import get_switchable_as_dense
from pypsa.opt import LConstraint, l_constraint, LExpression
from pyomo.environ import Var, NonNegativeReals


# adapted from pypsa.opf.extract_optimisation_results
def get_values(indexedvar):
    return pd.Series(indexedvar.get_values())


# adapted from pypsa.opf.extract_optimisation_results
def set_from_series(df, series, snapshots):
    df.loc[snapshots] = series.unstack(0).reindex(columns=df.columns)


def define_loss_constraints(network, snapshots, num_intervals=3):

    positions = range(1, num_intervals + 1)

    passive_branches = network.passive_branches()

    s_max_pus = get_switchable_as_dense(network, "Line", "s_max_pu")

    network.model.loss = Var(
        list(passive_branches.index), snapshots, domain=NonNegativeReals
    )

    loss_upper = {}
    loss_tangents = {}

    for branch in passive_branches.index:

        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        bt = branch[0]
        bn = branch[1]

        r_pu_eff = passive_branches.at[branch, "r_pu_eff"]

        s_nom_extendable = passive_branches.at[branch, "s_nom_extendable"]
        attr = "s_nom_max" if s_nom_extendable else "s_nom"
        s_nom_max = passive_branches.at[branch, attr]

        for sn in snapshots:

            s_max_pu = s_max_pus.loc[bn, sn]

            # adjust kcl
            # use of ._body because of pyomo bug
            for bus in [bus0, bus1]:
                network.model.power_balance[bus, sn]._body -= (
                    network.model.loss[bt, bn, sn] / 2
                )

            # adjust flow limits
            network.model.flow_upper[bt, bn, sn]._body += network.model.loss[bt, bn, sn]

            # upper loss limit
            lhs = LExpression(
                [(1, network.model.loss[bt, bn, sn])],
                -r_pu_eff * (s_max_pu * s_nom_max) ** 2,
            )
            loss_upper[bt, bn, sn] = LConstraint(lhs, "<=", LExpression())

            # loss tangents
            for k in positions:

                p_k = k / num_intervals * s_max_pu * s_nom_max
                loss_k = r_pu_eff * p_k ** 2
                slope_k = 2 * r_pu_eff * p_k
                offset_k = loss_k - slope_k * p_k

                for sign in (-1, 1):

                    lhs = LExpression([(1, network.model.loss[bt, bn, sn])])
                    rhs = LExpression(
                        [(sign * slope_k, network.model.passive_branch_p[bt, bn, sn])],
                        offset_k,
                    )
                    loss_tangents[sign, k, bt, bn, sn] = LConstraint(lhs, ">=", rhs)

    l_constraint(
        network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots
    )

    l_constraint(
        network.model,
        "loss_tangents",
        loss_tangents,
        list(positions),
        list(passive_branches.index),
        snapshots,
    )


def store_losses(network, snapshots, duals):

    network.lines_t["loss"] = pd.DataFrame(
        0, index=snapshots, columns=network.lines.index
    )

    loss_values = get_values(network.model.loss)

    set_from_series(network.lines_t.loss, loss_values.loc["Line"], snapshots)
