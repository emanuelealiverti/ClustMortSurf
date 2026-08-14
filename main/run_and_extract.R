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
sim_data <- simulate_scenario_4_1()
Y <- apply(sim_data$log_m, 1, function(x) t(x), simplify = FALSE)

# Model Fitting
out_MCMC <- run_model(
  Y = Y,
  n_iter = 2e4, print_step = 100,
  # Hyperparameters: alpha (tRPM persistence)
  a_alpha = 0.01, 
  b_alpha = 1,
  # Hyperparameters: delta_j^2 (GP prior variance)
  a_delta = 0.001, 
  b_delta = 0.001,
  # Hyperparameters: omega_j^2 (Innovation variance)
  a_omega = 0.001, 
  b_omega = 0.001,
  # Output path
  name_save = "res_simstudy_scenario1.RDS"
)



## ----------------------------------------------------------------------
## Extract results
## ----------------------------------------------------------------------
require(reshape2)
require(dplyr)

burnin = 1:1e4
true_clusters <- sim_data$clusters
p <- length(true_clusters)
TT <- ncol(true_clusters[[1]])

# Inference on partition
psm <- compute_prob_coclust(import = FALSE, save = FALSE, draws = out_MCMC$res$labels, burnin = burnin)
ppe <- compute_partition_point_est(psm = psm, draws = out_MCMC$res$labels, save = FALSE, burnin = burnin)
# Summarise partition
summary_partition <- summarise_clusters(ppe = ppe, psm = psm, clusters_true = clusters)
summary_partition$psm_acc
# Prepare Tex table (Table 1 in the manuscript)
toTex <- matrix(round(summary_partition$psm_acc, 3), nrow = p, ncol = TT)
rownames(toTex) <- paste0("Spline ", 1:p, " $(j = ", 1:p, ")$")
colnames(toTex) <- paste0("$t=", 1:TT, "$")
print(xtable::xtable(toTex, align = rep("c", 11)), sanitize.text.function = identity)


# Inference on individual coefficients
beta_units <- compute_beta_units(draws_beta = out_MCMC$res$beta, draws_labs = out_MCMC$res$labels, burnin = burnin)
beta_df <- summarise_beta(beta_units, sim_data$beta_units)
