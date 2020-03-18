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

import logging

logger = logging.getLogger(__name__)

# import loss models from ./loss_models.py
from _loss_models import *

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import *

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)


def remove_kvl(network, snapshots):

    formulation = snakemake.config["solving"]["options"].get("formulation", "kirchhoff")

    if formulation in ["angles", "cycles", "ptdf"]:
        n.model.del_component(network.model.passive_branch_p_def)
    
    if formulation in ["cycles", "kirchhoff"]:
        n.model.del_component(network.model.cycle_constraints)


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    config = snakemake.config
    model = snakemake.wildcards.model

    assert model in [
        "transport"
        "lossless",
        "lossy",
    ], f"The model {loss} has not been defined. Choose 'transport', 'lossless' or 'lossy'."

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:

        n = pypsa.Network(snakemake.input[0])

        n.lines.s_nom_max = n.lines.s_nom + config["additional_s_nom"]
        n.links.p_nom_max = config["links_p_nom_max"]

        n.lines = n.lines.loc[n.lines.s_nom != 0]

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])

        # set extra_functionality
        if model == "transport":
            extra_functionality = remove_kvl
        elif model == "lossy":
            extra_functionality = define_loss_constraints
        else:
            extra_functionality = None

        # set extra_postprocessing
        if model == "lossy":
            extra_postprocessing = store_losses
        else:
            extra_postprocessing = None


        n = solve_network(
            n,
            config=snakemake.config["solving"],
            solver_log=snakemake.log.solver,
            opts=snakemake.wildcards.opts,
            extra_functionality=extra_functionality,
            extra_postprocessing=extra_postprocessing
        )

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
