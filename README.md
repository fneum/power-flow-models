# Comparison of Power Flow and Loss Models in PyPSA-Eur

With rising shares of renewables and the need to properly
assess trade-offs between transmission, storage and sectoral integration as balancing options,
building a bridge between energy system models and detailed power flow
studies becomes increasingly important, but is computationally challenging.

In this paper, we compare both common and improved
approximations for two nonlinear phenomena,
power flow and transmission losses, in linear capacity expansion problems
that co-optimise investments in generation, storage and transmission infrastructure.

We evaluate different flow representations discussing differences in investment decisions,
nodal prices, the deviation of optimised flows and losses
from simulated AC power flows, and the computational performance.

By using the open European power system model \mbox{PyPSA-Eur}, 
that combines high spatial and temporal resolution,
we obtain detailed and reproducible benchmarks aiming at
facilitating the selection of a suitable power flow model.

Given the differences in complexity, the optimal choice
depends on the application, the access to computational
resources, and the level of spatial detail considered.

Although the commonly used transport model can already identify key features
of a cost-efficient system while being computationally performant,
deficiencies under high loading conditions are revealed
due to the lack of a physical grid representation.
Moreover, disregarding transmission losses overestimates optimal grid expansion by 20\%.

Adding a convex relaxation of quadratic losses with two or three tangents to the linearised
power flow equations and accounting for changing line impedances as the network is reinforced
suffices to accurately represent active power flows and losses in design studies.
These outputs are then sufficiently physical to be used in more detailed nonlinear simulations,
for instance to determine reactive power flows and voltages or dynamic analyses.

## Usage

Install common `pypsa-eur` environment.

```sh
conda activate pypsa-eur
cd power-flow-models
snakemake -j 99 check_all_powerflows
```