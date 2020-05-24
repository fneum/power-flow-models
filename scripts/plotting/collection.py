"""
Collection of plotting functions.
"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = "Copyright 2019-2020 Fabian Neumann (KIT), GNU GPL 3"

from .map import plot_network
from .stats import (
    check_capacities,
    check_costs,
    check_curtailment,
    check_energy_balance,
    check_energy_generated,
    check_energy_transmitted,
    check_flow_errors,
    check_slack,
)
from .single_lopf import plot_flow_vs_loss, plot_negative_marginal_prices
from .single_pf import plot_network_losses, plot_v_ang_diff
from .single_lopf_pf import (
    plot_duration_curve,
    plot_flow_comparison,
    plot_loss_comparison,
)
from .multiple_lopf import plot_performance, plot_cost_bar, plot_capacity_correlation, plot_price_duration_curve
from .convergence import plot_nonconverged, convergence_share
from .space import plot_feasible_space
