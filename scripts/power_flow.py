"""
Run power flow on solved networks.
"""

__author__ = "Fabian Neumann (KIT), Anika Bitsch (KIT)"
__copyright__ = (
    "Copyright 2019-2020 Fabian Neumann (KIT), Anika Bitsch (KIT), GNU GPL 3"
)

import pypsa
import pandas as pd

import logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    slack = True if snakemake.wildcards.slack == "distributed" else False

    n = pypsa.Network(snakemake.input[0])

    # remove load shedding generators (for distribution of slack)
    n.mremove("Generator", n.generators.loc[n.generators.carrier == "load"].index)

    set_components = n.controllable_one_port_components.union(
        n.controllable_branch_components
    ) - {"Load"}
    for c in n.iterate_components(set_components):
        if not c.df.empty:
            attr = "p0" if c.name == "Link" else "p"
            c.pnl.p_set = c.pnl.p_set.reindex(columns=c.df.index)
            c.pnl.p_set = c.pnl[attr]

    # set all buses to PV, since we don't know what Q set points are
    n.generators.control = "PV"

    # Need some PQ buses so that Jacobian doesn't break
    init_bus = n.buses.index[0]
    pq_gen_selection = n.generators[n.generators.bus == init_bus]
    n.generators.loc[pq_gen_selection.index, "control"] = "PQ"

    log = n.pf(distribute_slack=slack)

    pd.concat(log, axis=1).to_csv(snakemake.output.pf_log)

    n.export_to_netcdf(snakemake.output.network)
