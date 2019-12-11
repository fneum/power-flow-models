from os.path import normpath
configfile: "config.yaml"
subworkflow pypsaeur:
    workdir: "pypsa-eur"
    configfile: "config.pypsaeur.yaml"



COSTS="data/costs.csv"

wildcard_constraints:
    ll="(v|c)([0-9\.]+|opt|all)|all", # line limit, can be volume or cost
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    sectors="[+a-zA-Z0-9]+",
    opts="[-+a-zA-Z0-9\.]*",
    loss="[-+a-zA-Z0-9]+"

#rule solve_all_elec_lossy_networks:
    #input:pypsaeur(expand("results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",**config['scenario']))
    #input:expand("results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc",**config['scenario'])

def memory(w):
    factor = 10.
    #factor = 3.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))
        # return 4890+310 * int(w.clusters)

rule solve_lossy_network:
    input: pypsaeur("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc")
    output: 
      nc="results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc",
      csv=directory("results/csv/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}")
    log:
        solver="logs/{network}_s{simpl}_{clusters}_l{ll}_{opts}_{loss}_solver.log",
        python="logs/{network}_s{simpl}_{clusters}_l{ll}_{opts}_{loss}_python.log",
        memory="logs/{network}_s{simpl}_{clusters}_l{ll}_{opts}_{loss}_memory.log"
    threads: 4
    resources: mem=memory
    script: "scripts/solve_network_lossy.py"

rule lossy_power_flow_dist: 
    input: "results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc"
    output: nc ="results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}_dist_slack",
            csv=directory("results/csv/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}_dist_slack")
    script: "scripts/dist_slack.py"

rule lossy_power_flow: 
    input: "results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc"
    output: nc = "results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}_nr",
            csv=directory("results/csv/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}_nr")
    script: "scripts/newton_raphson.py"

rule power_flow_dist: 
    input: pypsaeur("results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc")
    output: nc ="results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_dist_slack",
            csv=directory("results/csv/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_dist_slack")
    script: "scripts/dist_slack.py"

rule power_flow: 
    input: pypsaeur("results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc")
    output: nc ="results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_nr",
            csv=directory("results/csv/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_nr")
    script: "scripts/newton_raphson.py"

rule plot_lossy_network:
    input:
        network="results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc",
        tech_costs=COSTS
    output:
        only_map="results/plots/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{loss}.{ext}",
        ext="results/plots/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext_{loss}.{ext}"
    script: "scripts/plot_network.py"

def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return ([COSTS] +
            expand("results/networks/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{loss}.nc",
                   network=w.network,
                   ll=ll,
                   **{k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
                      for k in ["simpl", "clusters", "opts"]}))

rule make_summary:
    input: input_make_summary
    output: directory("results/summaries/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{loss}")
    script: "scripts/make_summary.py"

rule plot_summary:
    input: "results/summaries/{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{loss}"
    output: "results/plots/summary_{summary}_{network}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{loss}.{ext}"
    script: "scripts/plot_summary.py"

def input_plot_p_nom_max(wildcards):
    return [('networks/{network}_s{simpl}{maybe_cluster}_{loss}.nc'
             .format(maybe_cluster=('' if c == 'full' else ('_' + c)), **wildcards))
            for c in wildcards.clusters.split(",")]
rule plot_p_nom_max:
    input: input_plot_p_nom_max
    output: "results/plots/{network}_s{simpl}_cum_p_nom_max_{clusters}_{technology}_{country}_{loss}.{ext}"
    script: "scripts/plot_p_nom_max.py"
