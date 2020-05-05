"""
Plotting functions regarding the NR power flow convergence rate of subnetworks.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

import matplotlib.pyplot as plt


# careful, this is a bit hard-coded!
def convergence_share(log, fm):
    conv = log.converged.astype(int)
    conv.columns = [
        "Continental Europe",
        "Nordic",
        "Baltic",
        "Mallorca",
        "Ireland",
        "Great Britain",
        "Sicily",
    ]
    conv = conv.drop(columns=["Mallorca", "Sicily"])
    conv_share = (1 - conv.sum() / len(conv)) * 100  # %
    conv_share.name = fm
    return conv_share


def plot_nonconverged(conv_share, model_names, fn=None):

    fig, ax = plt.subplots(figsize=(6, 2.5))
    conv_share.rename(columns=model_names, inplace=True)
    conv_share.T.plot.bar(ax=ax)

    plt.legend(title="Synchronous Zone")
    plt.ylabel("Snapshots not\nconverged [%]")

    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")
