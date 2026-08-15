#########################################################################
## Execution script for Simulation Study in Section 4.1                ##
## of "Bayesian local clustering of age-period mortality surfaces      ##
## across multiple countries" (Romanò, Aliverti and Durante)           ##
#########################################################################

## ----------------------------------------------------------------------
## Run Gibbs sampler
## ----------------------------------------------------------------------

rm(list = ls())

# Load Utilities & MCMC Routine
source("../main/simulate_data.R")
source("../main/utils.R")
source("../main/gibbs.R")

# Data Preparation
n  <- 5 # number of countries
TT <- 10 # number of time points
p  <- 6 # number of spline bases

sim      <- simulate_scenario_4_1(n = n, TT = TT, p = p)
clusters <- sim$clusters

# Model Fitting
out_MCMC <- run_model(Y = sim$log_m, ages = sim$ages,
  n_iter = 2e4, print_step = 500,
  a_alpha = 0.01, b_alpha = 1, a_delta = 0.001, b_delta = 0.001, 
  a_omega = 0.001, b_omega = 0.001,
  name_save = "res_simstudy_scenario1.RDS")



## ----------------------------------------------------------------------
## Extract results
## ----------------------------------------------------------------------
burnin = 1:1e4

# Inference on partition
psm <- compute_prob_coclust(import = FALSE, save = FALSE, draws = out_MCMC$res$labels, burnin = burnin)
ppe <- compute_partition_point_est(psm = psm, draws = out_MCMC$res$labels, save = FALSE, burnin = burnin)
coclust <- list(psm = psm, ppe = ppe, true = clusters)

# Inference on individual coefficients
beta_est <- compute_beta_units(draws_beta = out_MCMC$res$beta, draws_labs = out_MCMC$res$labels, burnin = burnin)
beta_countries <- summarise_beta(beta_est, sim$beta)

saveRDS(coclust, "coclust.RDS")
saveRDS(beta_countries, "beta_countries.RDS")