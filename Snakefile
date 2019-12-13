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
    output: 
      nc="results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc",
      csv=directory("results/csv/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}")
    log:
        solver="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_solver.log",
        python="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_python.log",
        memory="logs/elec_s_{clusters}_lcopt_{opts}_L{loss}_memory.log"
    threads: 4
    resources: mem=memory
    script: "scripts/solve_network_lossy.py"

rule solve_all_lossy_networks:
    input: 
        expand("results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc",
               **config["scenario"])


# PF CHECK RULES

rule pf_for_lossy_network: 
    input: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc"
    output: nc="results/pf/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}_S{slack}.nc",
            csv=directory("results/csv/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}_S{slack}")
    script: "scripts/check_power_flow.py"

rule pf_for_lossless_network: 
    input: pypsaeur("results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: nc="results/pf/elec_s_{clusters}_ec_lcopt_{opts}_S{slack}.nc",
            csv=directory("results/csv/elec_s_{clusters}_ec_lcopt_{opts}_S{slack}")
    script: "scripts/check_power_flow.py"

rule pf_for_all:
    input: 
        lossy=expand("results/pf/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}_S{slack}.nc",
                     **config["scenario"]),
        lossless=expand("results/pf/elec_s_{clusters}_ec_lcopt_{opts}_S{slack}.nc",
                     **config["scenario"])

# EVALUATION

# COSTS="data/costs.csv"

# rule plot_lossy_network:
#     input:
#         network="results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc",
#         tech_costs=COSTS
#     output:
#         only_map="results/plots/elec_s_{clusters}_ec_lcopt_{opts}_{attr}_L{loss}.{ext}",
#         ext="results/plots/elec_s_{clusters}_ec_lcopt_{opts}_{attr}_ext_L{loss}.{ext}"
#     script: "scripts/plot_network.py"

# def input_make_summary(w):
#     # It's mildly hacky to include the separate costs input as first entry
#     if w.ll.endswith("all"):
#         ll = config["scenario"]["ll"]
#         if len(w.ll) == 4:
#             ll = [l for l in ll if l[0] == w.ll[0]]
#     else:
#         ll = w.ll
#     return ([COSTS] +
#             expand("results/networks/elec_s_{clusters}_ec_lcopt_{opts}_L{loss}.nc",
#                    network=w.network,
#                    ll=ll,
#                    **{k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
#                       for k in ["simpl", "clusters", "opts"]}))

# rule make_summary:
#     input: input_make_summary
#     output: directory("results/summaries/elec_s_{clusters}_ec_lcopt_{opts}_{country}_L{loss}")
#     script: "scripts/make_summary.py"

# rule plot_summary:
#     input: "results/summaries/elec_s_{clusters}_ec_lcopt_{opts}_{country}_L{loss}"
#     output: "results/plots/summary_{summary}_elec_s_{clusters}_ec_lcopt_{opts}_{country}_L{loss}.{ext}"
#     script: "scripts/plot_summary.py"

# def input_plot_p_nom_max(wildcards):
#     return [('networks/elec_s{maybe_cluster}_L{loss}.nc'
#              .format(maybe_cluster=('' if c == 'full' else ('_' + c)), **wildcards))
#             for c in wildcards.clusters.split(",")]
            
# rule plot_p_nom_max:
#     input: input_plot_p_nom_max
#     output: "results/plots/elec_s_cum_p_nom_max_{clusters}_{technology}_{country}_L{loss}.{ext}"
#     script: "scripts/plot_p_nom_max.py"
