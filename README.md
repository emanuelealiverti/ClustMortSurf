# Bayesian local clustering of age-period mortality surfaces across multiple countries

This repository contains the code implementing the model, along with the materials to reproduce the
simulation study, of the paper

> Romanò, G., Aliverti, E. and Durante, D. (2025+). *Bayesian local clustering of age-period
> mortality surfaces across multiple countries.*
> [arXiv:2504.05240](https://arxiv.org/abs/2504.05240)

The paper develops a Bayesian model for multi-country log-mortality rates in which the age structure
of mortality is represented through a **B-spline expansion**, whose country-specific dynamic
coefficients are endowed with a **temporal random partition prior** (tRPM; Page et al., 2022). This
construction lets countries cluster **locally**, i.e. differently across combinations of age classes
and periods, rather than globally, and learns the number of groups automatically.

## Repository structure

```
.
├── main/                    core functions
│   ├── simulate_data.R      B-spline bases and data-generating mechanism
│   └── gibbs_tRPM.R         (to be added) Gibbs sampler and Monte Carlo summaries
└── Simulations/             reproducible tutorial for the simulation study
    ├── Simulation_4-1.md    step-by-step tutorial (Section 4.1)
    ├── beta_df.RDS          posterior means of the spline coefficients
    ├── coclust_df.RDS       (to be added) estimated cluster memberships
    └── img/                 figures produced by the tutorial
```

- **[`main/`](main)** collects the core functions: the Gibbs sampling algorithm derived in Section 3.1 of the manuscript and the data-generating mechanism of the simulation study 
- **[`Simulations/`](Simulations)** contains the tutorial
  [`Simulation_4-1.md`](Simulations/Simulation_4-1.md), reproducing the Scenario proposed in Section 4.1 of the manuscript


## References

- Page, G. L., Quintana, F. A. and Dahl, D. B. (2022). Dependent modeling of temporal sequences of random partitions. *Journal of Computational and Graphical Statistics*, 31, 614–627.
