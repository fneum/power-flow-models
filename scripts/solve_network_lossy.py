import pypsa
import numpy as np
import pandas as pd
import os
import sys
import gc
import logging
import pyomo.kernel as pmo 
import pyomo.environ as penv

from pypsa import opf
from pypsa.opt import LConstraint, l_constraint, LExpression
from pypsa.descriptors import free_output_series_dataframes
from six import iteritems, itervalues, string_types
from pyomo.environ import Constraint, Objective, Var, ComponentUID
from vresutils.benchmark import memory_logger

import logging
logger = logging.getLogger(__name__)

# import loss models from ./loss_models.py
from _loss_models import *

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd()+"/pypsa-eur/scripts")

from solve_network import *

# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)


if __name__ == "__main__":

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        patch_pyomo_tmpdir(tmpdir)

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    config = snakemake.config

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        
        n = pypsa.Network(snakemake.input[0])

        n.lines.s_nom_max = n.lines.s_nom + config["additional_s_nom"]
        n.links.p_nom_max = config["links_p_nom_max"]

        n = prepare_network(n)
        
        try:
            n = solve_network(n, opts=snakemake.config['solving'],
                              solver_logs=snakemak.logs.solver,
                              extra_functionality=globals()[f"{loss}"],
                              extra_postprocessing=post_processing)
        except KeyError:
            raise RuntimeError(f"The loss function {loss} has not been defined.")
            
        n.export_to_netcdf(snakemake.output.nc)
        n.export_to_csv_folder(snakemake.output.csv)

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))