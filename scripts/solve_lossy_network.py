"""
Solve lossy networks.
"""

__author__ = "Fabian Neumann (KIT), Anika Bitsch (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), Anika Bitsch (KIT), GNU GPL 3"

import pypsa
import numpy as np
import pandas as pd
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


if __name__ == "__main__":

    tmpdir = snakemake.config["solving"].get("tmpdir")
    if tmpdir is not None:
        patch_pyomo_tmpdir(tmpdir)

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    config = snakemake.config
    loss = snakemake.wildcards.loss

    assert loss in [
        "cosine",
        "square",
    ], f"The loss function {loss} has not been defined. Try 'cosine' or 'square'"

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:

        n = pypsa.Network(snakemake.input[0])

        n.lines.s_nom_max = n.lines.s_nom + config["additional_s_nom"]
        n.links.p_nom_max = config["links_p_nom_max"]

        n.lines = n.lines.loc[n.lines.s_nom != 0]

        n = prepare_network(n, solve_opts=snakemake.config["solving"]["options"])

        n = solve_network(
            n,
            config=snakemake.config["solving"],
            solver_log=snakemake.log.solver,
            opts=snakemake.wildcards.opts,
            extra_functionality=globals()[f"{loss}"],
            extra_postprocessing=post_processing,
        )

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
