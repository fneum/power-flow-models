configfile: "config.yaml"

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    configfile: "config.pypsaeur.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    opts="[-+a-zA-Z0-9\.]*",
    model="(transport|lossless-[0-9]+|lossy-[0-9]+-[0-9]+)",
    slack="(distributed|regular)"


# adapted from https://github.com/PyPSA/pypsa-eur/blob/master/Snakefile
def memory(w):
    factor = 3.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    if "lossy" in w.model:
        factor *= 0.33 * int(w.model.split("-")[-1])
    return int(factor * (10000 + 195 * int(w.clusters)))


# SOLVING RULES

rule solve_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_M{model}.nc"
    log:
        solver="logs/elec_s_{clusters}_lcopt_{opts}_M{model}_solver.log",
        python="logs/elec_s_{clusters}_lcopt_{opts}_M{model}_python.log",
        memory="logs/elec_s_{clusters}_lcopt_{opts}_M{model}_memory.log"
    threads: 4
    resources: mem=memory
    script: "scripts/solve_network.py"

rule solve_all_networks:
    input: 
        expand("results/networks/elec_s_{clusters}_ec_lcopt_{opts}_M{model}.nc",
               **config["scenario"])


# POWER FLOW RULES

rule check_powerflow: 
    input: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_M{model}.nc"
    log:
        python="logs/elec_s_{clusters}_lcopt_{opts}_M{model}_S{slack}_python.log"
    output:
        network="results/pf/elec_s_{clusters}_ec_lcopt_{opts}_M{model}_S{slack}.nc",
        pf_log="results/pf/log_elec_s_{clusters}_ec_lcopt_{opts}_M{model}_S{slack}.csv"
    script: "scripts/power_flow.py"

rule check_all_powerflows:
    input: 
        expand("results/pf/elec_s_{clusters}_ec_lcopt_{opts}_M{model}_S{slack}.nc",
               **config["scenario"]),
