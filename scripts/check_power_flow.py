import pypsa

import logging
logger = logging.getLogger(__name__)


if __name__ == "__main__":

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    slack = True if snakemake.wildcards.slack == 'distributed' else False

    n = pypsa.Network(snakemake.input[0])

    n.generators_t.p_set = n.generators_t.p_set.reindex(columns=n.generators.index)
    n.generators_t.p_set = n.generators_t.p
    n.storage_units_t.p_set = n.storage_units_t.p_set.reindex(columns=n.storage_units.index)
    n.storage_units_t.p_set = n.storage_units_t.p

    #set all buses to PV, since we don't know what Q set points are
    n.generators.control = "PV"

    #Need some PQ buses so that Jacobian doesn't break
    init_bus = n.buses.index[0]
    pq_gen_selection = n.generators[n.generators.bus==init_bus] # TODO select first bus
    n.generators.loc[pq_gen_selection.index,"control"] = "PQ"
    
    n.pf(distribute_slack=slack)
    
    n.export_to_netcdf(snakemake.output.nc)
    n.export_to_csv_folder(snakemake.output.csv)
