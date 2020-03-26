"""
Extra functionality implementing loss models.
"""

__author__ = "Fabian Neumann (KIT), Anika Bitsch (KIT)"
__copyright__ = (
    "Copyright 2019-2020 Fabian Neumann (KIT), Anika Bitsch (KIT), GNU GPL 3"
)

import pandas as pd
import numpy as np

from pypsa.descriptors import get_switchable_as_dense
from pypsa.opt import LConstraint, l_constraint, LExpression
from pyomo.environ import Var, NonNegativeReals


# adapted from pypsa.opf.define_passive_branch_constraints
def redo_passive_branch_constraints(network, snapshots):

    model_components_to_delete = [
        "flow_upper",
        "flow_lower",
        "flow_upper_index",
        "flow_lower_index",
        "flow_upper_index_0",
        "flow_lower_index_0",
        "flow_upper_index_1",
        "flow_lower_index_1",
    ]
    for model_component in model_components_to_delete:
        network.model.del_component(model_component)

    passive_branches = network.passive_branches()
    extendable_branches = passive_branches[passive_branches.s_nom_extendable]
    fixed_branches = passive_branches[~passive_branches.s_nom_extendable]

    s_max_pu = pd.concat(
        {
            c: get_switchable_as_dense(network, c, "s_max_pu", snapshots)
            for c in network.passive_branch_components
        },
        axis=1,
        sort=False,
    )

    flow_upper = {
        (b[0], b[1], sn): [
            [
                (1, network.model.passive_branch_p[b[0], b[1], sn]),
                (1, network.model.loss[b[0], b[1], sn]),
            ],
            "<=",
            s_max_pu.at[sn, b] * fixed_branches.at[b, "s_nom"],
        ]
        for b in fixed_branches.index
        for sn in snapshots
    }

    flow_upper.update(
        {
            (b[0], b[1], sn): [
                [
                    (1, network.model.passive_branch_p[b[0], b[1], sn]),
                    (1, network.model.loss[b[0], b[1], sn]),
                    (
                        -s_max_pu.at[sn, b],
                        network.model.passive_branch_s_nom[b[0], b[1]],
                    ),
                ],
                "<=",
                0,
            ]
            for b in extendable_branches.index
            for sn in snapshots
        }
    )

    l_constraint(
        network.model, "flow_upper", flow_upper, list(passive_branches.index), snapshots
    )

    flow_lower = {
        (b[0], b[1], sn): [
            [
                (1, network.model.passive_branch_p[b[0], b[1], sn]),
                (-1, network.model.loss[b[0], b[1], sn]),
            ],
            ">=",
            -s_max_pu.at[sn, b] * fixed_branches.at[b, "s_nom"],
        ]
        for b in fixed_branches.index
        for sn in snapshots
    }

    flow_lower.update(
        {
            (b[0], b[1], sn): [
                [
                    (1, network.model.passive_branch_p[b[0], b[1], sn]),
                    (-1, network.model.loss[b[0], b[1], sn]),
                    (
                        s_max_pu.at[sn, b],
                        network.model.passive_branch_s_nom[b[0], b[1]],
                    ),
                ],
                ">=",
                0,
            ]
            for b in extendable_branches.index
            for sn in snapshots
        }
    )

    l_constraint(
        network.model, "flow_lower", flow_lower, list(passive_branches.index), snapshots
    )


# adapted from pypsa.opf.extract_optimisation_results
def get_values(indexedvar):
    return pd.Series(indexedvar.get_values())


# adapted from pypsa.opf.extract_optimisation_results
def set_from_series(df, series, snapshots):
    df.loc[snapshots] = series.unstack(0).reindex(columns=df.columns)


def define_loss_constraints(network, snapshots):

    tangents = network.tangents

    positions = range(1, tangents + 1)
    signs = [-1, 1]

    passive_branches = network.passive_branches()

    s_max_pus = get_switchable_as_dense(network, "Line", "s_max_pu")

    network.model.loss = Var(
        list(passive_branches.index), snapshots, domain=NonNegativeReals
    )

    redo_passive_branch_constraints(network, snapshots)

    loss_upper = {}
    loss_tangents = {}

    for branch in passive_branches.index:

        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        bt = branch[0]
        bn = branch[1]

        r_pu_eff = passive_branches.at[branch, "r_pu_eff"]

        if passive_branches.at[branch, "s_nom_extendable"]:
            attr = "s_nom_max"
        elif passive_branches.at[branch, "s_nom_opt"] != 0.0:
            attr = "s_nom_opt"
        else:
            attr = "s_nom"

        s_nom_max = passive_branches.at[branch, attr]

        assert np.isfinite(s_nom_max) and not np.isnan(
            s_nom_max
        ), f"Infinite or non-existent 's_nom_max' encountered at line {bn}"

        for sn in snapshots:

            s_max_pu = s_max_pus.loc[sn, bn]

            # adjust kcl
            # use of ._body because of pyomo bug
            for bus in [bus0, bus1]:
                network.model.power_balance[bus, sn]._body -= (
                    network.model.loss[bt, bn, sn] / 2
                )

            # upper loss limit
            lhs = LExpression(
                [(1, network.model.loss[bt, bn, sn])],
                -r_pu_eff * (s_max_pu * s_nom_max) ** 2,
            )
            loss_upper[bt, bn, sn] = LConstraint(lhs, "<=", LExpression())

            # loss tangents
            for k in positions:

                p_k = k / tangents * s_max_pu * s_nom_max
                loss_k = r_pu_eff * p_k ** 2
                slope_k = 2 * r_pu_eff * p_k
                offset_k = loss_k - slope_k * p_k

                for sign in signs:

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
        signs,
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
