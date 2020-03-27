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

from vresutils.benchmark import memory_logger
from pypsa.opt import l_constraint

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
    
    if not hasattr(n.links, 'reversed'): return
    
    ext_rev_links = network.links.loc[
                        (network.links.reversed==True) & 
                        (network.links.p_nom_extendable==True)
                    ].index

    if len(ext_rev_links) == 0: return
    
    constraints = {lk :
            [[(1,network.model.link_p_nom[lk.split("-")[0]]),
              (-1,network.model.link_p_nom[lk])
             ],"==",0.]
        for lk in ext_rev_links}
    
    l_constraint(network.model, "bidirectional_link", constraints, list(ext_rev_links))


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    config = snakemake.config

    model_wc = snakemake.wildcards.model.split("-")
    model = model_wc[0]

    assert model in [
        "transport",
        "lossless",
        "lossy",
    ], f"The model {model} has not been defined. Choose 'transport', 'lossless' or 'lossy'."

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:

        n = pypsa.Network(snakemake.input[0])

        lk_config = config["links"]
        n.lines.s_nom_max = n.lines.apply(
            lambda line: max(
                line.s_nom + lk_config["s_nom_add"],
                line.s_nom * lk_config["s_nom_factor"],
            ),
            axis=1
        )
        n.lines = n.lines.loc[n.lines.s_nom != 0]

        lk_config = config["lines"]
        n.links.p_nom_max = lk_config["p_nom_max"]
        n.links.efficiency = n.links.apply(lambda lk: 1 - lk.length * lk_config["efficiency_per_length"], axis=1)
        split_bidirectional_links(n)

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])

        # set iterating
        if model == "transport":
            skip_iterating = True
        elif model == "lossy":
            n.tangents = int(model_wc[1])
            skip_iterating = False
        elif model == "lossless":
            iterations = int(model_wc[1])
            if iterations == 0:
                skip_iterating = True
            else:
                skip_iterating = False
                snakemake.config["solving"]["options"]["min_iterations"] = iterations - 1
                snakemake.config["solving"]["options"]["max_iterations"] = iterations

        def extra_functionality(network, snapshots):
            tie_bidirectional_link_p_nom(network, snapshots)
            if model == "transport":
                remove_kvl_constraints(netowrk, snapshots)
            if model == "lossy"
                define_loss_constraints(network, snapshots)

        def extra_postprocessing(network, snapshots, duals):
            if model == "lossy":
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

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
