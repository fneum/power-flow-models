

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake, Dict
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='45', lv='1.0', opts='Co2L-3H'),
            input=["networks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc"],
            output=["results/networks/s{simpl}_{clusters}_lv{lv}_{opts}.nc"],
            log=dict(solver="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_solver.log",
                     python="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_python.log")
        )

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        patch_pyomo_tmpdir(tmpdir)

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])
        n.generators_t.p_set = n.generators_t.p_set.reindex(columns=n.generators.index)
        n.generators_t.p_set = n.generators_t.p
        n.storage_units_t.p_set = n.storage_units_t.p_set.reindex(columns=n.storage_units.index)
        n.storage_units_t.p_set = n.storage_units_t.p

        #set all buses to PV, since we don't know what Q set points are
        n.generators.control = "PV"

        #Need some PQ buses so that Jacobian doesn't break
        f = n.generators[network.generators.bus == "492"]
        n.generators.loc[f.index,"control"] = "PQ"
        n.pf()
        n.export_to_netcdf(snakemake.output.nc)
        n.export_to_csv_folder(snakemake.output.csv)

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))