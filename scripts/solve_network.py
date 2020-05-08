"""
Solve lossy networks.
"""

__author__ = "Fabian Neumann (KIT), Anika Bitsch (KIT)"
__copyright__ = (
    "Copyright 2019-2020 Fabian Neumann (KIT), Anika Bitsch (KIT), GNU GPL 3"
)

import pypsa
import os
import sys
import numpy as np
import pandas as pd

from vresutils.benchmark import memory_logger
from pypsa.opt import l_constraint
from copy import deepcopy

import logging

logger = logging.getLogger(__name__)

# import loss models from ./loss_models.py
from _loss_models import *

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import *

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)


def split_bidirectional_links(n):

    n.links.p_min_pu = 0
    n.links["reversed"] = False
    rev_links = n.links.copy().rename({"bus0": "bus1", "bus1": "bus0"}, axis=1)
    rev_links.capital_cost = 0
    rev_links["reversed"] = True
    rev_links.index = [f"{i}-reversed" for i in rev_links.index]

    n.links = pd.concat([n.links, rev_links], sort=False)


def remove_kvl_constraints(network, snapshots):

    formulation = snakemake.config["solving"]["options"].get("formulation", "kirchhoff")

    if formulation in ["angles", "cycles", "ptdf"]:
        n.model.del_component(network.model.passive_branch_p_def)

    if formulation in ["cycles", "kirchhoff"]:
        n.model.del_component(network.model.cycle_constraints)


def tie_bidirectional_link_p_nom(network, snapshots):

    if not hasattr(n.links, "reversed"):
        return

    ext_rev_links = network.links.loc[
        (network.links.reversed == True) & (network.links.p_nom_extendable == True)
    ].index

    if len(ext_rev_links) == 0:
        return

    constraints = {
        lk: [
            [
                (1, network.model.link_p_nom[lk.split("-")[0]]),
                (-1, network.model.link_p_nom[lk]),
            ],
            "==",
            0.0,
        ]
        for lk in ext_rev_links
    }

    l_constraint(network.model, "bidirectional_link", constraints, list(ext_rev_links))

    network.model.bidirectional_link.pprint()


def update_line_parameters(n):

    lines_ext_b = n.lines.s_nom_extendable

    if not lines_ext_b.any():
        return

    lines = pd.DataFrame(n.lines[["r", "x", "type", "num_parallel"]])

    lines["s_nom"] = (
        np.sqrt(3)
        * n.lines["type"].map(n.line_types.i_nom)
        * n.lines.bus0.map(n.buses.v_nom)
    ).where(n.lines.type != "", n.lines["s_nom"])

    lines_ext_untyped_b = (n.lines.type == "") & lines_ext_b
    lines_ext_typed_b = (n.lines.type != "") & lines_ext_b

    if lines_ext_untyped_b.any():
        for attr in ("r", "x"):
            n.lines.loc[lines_ext_untyped_b, attr] = lines[attr].multiply(
                lines["s_nom"] / n.lines["s_nom_opt"]
            )

    if lines_ext_typed_b.any():
        n.lines.loc[lines_ext_typed_b, "num_parallel"] = (
            n.lines["s_nom_opt"] / lines["s_nom"]
        )


def modify_line_attributes(n, ln_config):
    n.lines.s_nom_max = n.lines.apply(
        lambda line: max(
            line.s_nom + ln_config["s_nom_add"], line.s_nom * ln_config["s_nom_factor"],
        ),
        axis=1,
    )
    n.lines = n.lines.loc[n.lines.s_nom != 0]
    n.lines.s_max_pu = ln_config["s_max_pu"]


def modify_link_attributes(n, lk_config):
    n.links.p_nom_max = lk_config["p_nom_max"]
    if "lossy" in n.flow_model:
        n.links.efficiency = n.links.apply(
            lambda lk: 1 - lk.length * lk_config["loss_per_length"], axis=1
        )


def ac_lines_to_links(n):

    lk_attrs = [
        "bus0",
        "bus1",
        "length",
        "s_nom",
        "s_max_pu",
        "capital_cost",
        "s_nom_max",
        "s_nom_min",
        "s_nom_extendable",
    ]
    lines_lk = n.lines.loc[:, lk_attrs]

    # translate s_ to p_
    lines_lk.columns = [
        f"p_{col[2:]}" if col.startswith("s_") else col for col in lines_lk.columns
    ]

    # need to ensure unique naming
    lines_lk.index = [f"AC_{ln}" for ln in lines_lk.index]

    lines_lk["efficiency"] = lines_lk.apply(lambda ln: 1 - ln.length * 3e-5, axis=1)

    lines_lk["p_min_pu"] = -lines_lk.p_max_pu
    lines_lk["carrier"] = "AC"

    n.import_components_from_dataframe(lines_lk, "Link")

    n.mremove("Line", n.lines.index)


def ac_links_to_lines(n, lines_orig):

    def _get_ac_links(n, reverse=True):

        links = n.links.copy()
        links_t_p0 = n.links_t.p0.copy()
        links_t_p1 = n.links_t.p1.copy()

        if reverse:
            ln = links.loc[(links.reversed == True) & (links.carrier == "AC")]
        else:
            ln = links.loc[(links.reversed != True) & (links.carrier == "AC")]

        ln_p0 = links_t_p0.loc[:, ln.index]
        ln_p1 = links_t_p1.loc[:, ln.index]

        def orig_name(l):
            return l.split("-")[0].split("_")[1]

        ln.index = [orig_name(l) for l in ln.index]
        ln_p0.columns = [orig_name(l) for l in ln_p0.columns]
        ln_p1.columns = [orig_name(l) for l in ln_p1.columns]

        return ln, ln_p0, ln_p1

    lines, p0, p1 = _get_ac_links(n, reverse=False)
    lines_rev, p0_rev, p1_rev = _get_ac_links(n, reverse=True)

    n.import_components_from_dataframe(lines_orig, "Line")

    n.lines.s_nom_opt = lines.p_nom_opt
    n.lines_t["p0"] = p1_rev - p1
    n.lines_t["p1"] = -n.lines_t.p0  # p0 and p1 store effective flow (received flow)
    n.lines_t["loss"] = p0 + p1 + p0_rev + p1_rev

    n.mremove("Link", n.links.loc[n.links.carrier == "AC"].index)
    n.links.drop(columns=["underground", "under_construction"], inplace=True) # these cause trouble when exporting as .nc

    n.calculate_dependent_values()
    n.determine_network_topology()


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    config = snakemake.config

    flow_model_wc = snakemake.wildcards.model.split("-")
    flow_model = flow_model_wc[0]

    assert flow_model in [
        "transport",
        "lossytransport",
        "lossless",
        "lossy",
    ], f"The flow model {flow_model} has not been defined. Choose 'transport', 'lossytransport', 'lossless' or 'lossy'."

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:

        n = pypsa.Network(snakemake.input[0])

        n.flow_model = flow_model

        modify_line_attributes(n, config["lines"])
        modify_link_attributes(n, config["links"])

        lines_orig = n.lines.copy()

        if flow_model == "lossytransport":
            ac_lines_to_links(n)

        split_bidirectional_links(n)

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])

        # set iterations
        if flow_model in ["transport", "lossytransport"]:
            skip_iterating = True
        elif flow_model in ["lossless", "lossy"]:
            iterations = int(flow_model_wc[-1])
            if iterations == 0:
                skip_iterating = True
            else:
                skip_iterating = False
                snakemake.config["solving"]["options"]["min_iterations"] = iterations
                snakemake.config["solving"]["options"]["max_iterations"] = iterations

        # set tangents
        if flow_model == "lossy":
            n.tangents = int(flow_model_wc[-2])

        def extra_functionality(network, snapshots):
            tie_bidirectional_link_p_nom(network, snapshots)
            if network.flow_model == "transport":
                remove_kvl_constraints(network, snapshots)
            elif network.flow_model == "lossy":
                define_loss_constraints(network, snapshots)

        def extra_postprocessing(network, snapshots, duals):
            if network.flow_model == "lossy":
                store_losses(network, snapshots, duals)

        n = solve_network(
            n,
            config=snakemake.config["solving"],
            solver_log=snakemake.log.solver,
            opts=snakemake.wildcards.opts,
            extra_functionality=extra_functionality,
            extra_postprocessing=extra_postprocessing,
            skip_iterating=skip_iterating,
        )

        if flow_model == "lossytransport":
            ac_links_to_lines(n, lines_orig)

        update_line_parameters(n)

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
