# Bayesian local clustering of age-period mortality surfaces across multiple countries

This repository contains the code implementing the model, along with the materials to reproduce the simulation study described in Section 4.1 of

> Romanò, G., Aliverti, E. and Durante, D. (2025+). *Bayesian local clustering of age-period mortality surfaces across multiple countries.* [arXiv:2504.05240](https://arxiv.org/abs/2504.05240)

The paper develops a Bayesian model for multi-country log-mortality surfaces, in which the age structure of mortality is represented through a B-spline expansion whose country-specific dynamic coefficients are endowed with a temporal random partition process. This construction allows countries cluster *locally*, i.e. differently across combinations of age classes and calendar-years.

## Repository structure

```
.
├── main/                    
│   ├── gibbs.R              Gibbs sampler functions
│   ├── run_and_extract.R    Run Gibbs sampler save and Monte Carlo summaries
│   └── utils.R              Utility functions for both Gibbs sampler and results extraction
└── Simulations/             
    ├── Simulation_4-1.md    step-by-step tutorial (Section 4.1)
    ├── simulate_data.R      Simulate articial log-mortality surfaces
    ├── beta_countries.RDS   posterior means of the spline coefficients (MCMC output)
    ├── coclust.RDS          posterior similarity matrix and estimated partitions (MCMC output)
    └── img/                 figures 
```

- **[`main/`](main)** collects the core functions implementing the Gibbs sampling algorithm (presented in Section 3.1 of the manuscript) and various utilities for post-processing
- **[`Simulations/`](Simulations)** contains the tutorial [`Simulation_4-1.md`](Simulations/Simulation_4-1.md), and the functions the MCMC output and the figures associated to the simulation study

