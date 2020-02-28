configfile: "config.yaml"

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    configfile: "config.pypsaeur.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9\.]*",
    loss="[-+a-zA-Z0-9]+",
    slack="(distributed|regular)"


def memory(w):
    factor = 6.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


# LOPF WITH LOSSY FORMULATIONS RULES

rule solve_lossy_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc"
    log:
        solver="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_solver.log",
        python="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_python.log",
        memory="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_memory.log"
    threads: 4
    resources: mem=memory
    script: "scripts/solve_lossy_network.py"

rule solve_all_lossy_networks:
    input: 
        expand("results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc",
               **config["scenario"])


# PF CHECK RULES

rule lossy_network_pf: 
    input: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc"
    output: "results/pf/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}_S{slack}.nc"
    script: "scripts/power_flow.py"

rule lossless_network_pf: 
    input: pypsaeur("results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/pf/elec_s_{clusters}_ec_lcopt_{opts}_S{slack}.nc"
    script: "scripts/power_flow.py"

rule all_pf:
    input: 
        lossy=expand("results/pf/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}_S{slack}.nc",
                     **config["scenario"]),
        lossless=expand("results/pf/elec_s_{clusters}_ec_lcopt_{opts}_S{slack}.nc",
                     **config["scenario"])
