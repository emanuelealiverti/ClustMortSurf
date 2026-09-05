

# Simulation study

This tutorial reproduces the simulation scenario presented in Section 4.1 of the manuscript [*Bayesian local clustering of age-period mortality surfaces across multiple countries*](https://arxiv.org/abs/2504.05240). The goal of the study is to generate mortality surfaces with complex dynamic clustering patterns, and assess to what extent the proposed formulation is able to recover (i) such complex local grouping structures among countries that vary across ages and periods, and (ii) the associated cluster-specific coefficients characterizing mortality rates.

We focus on artificial mortality rates for `n = 5` synthetic countries, observed over `T = 10` periods and for the ages `X = {0, 1, ..., 100}`.

The code is meant to be run from within this folder (\[`Simulations/`\]), and the core routines are sourced from [`main/`](../main).

``` r
# Gibbs sampler utilities
source("../main/gibbs.R")
source("../main/utils.R")
```

------------------------------------------------------------------------

## 1. Data-generating mechanism

Data are generated directly from the model proposed in the paper — as implemented in [`simulate_data.R`](simulate_data.R). Under the proposed specification, the age structure of mortality in each country and calendar year is constructed by a B-spline expansion, whose country-specific coefficients are **tied** whenever two countries belong to the same cluster at a given pair (basis `j`, calendar year `t`). In this simulated scenario, the local partitions are fixed at prespecified levels, centering the cluster-specific coefficients around simple linear trends in time.

The cluster memberships $c_{jt}$ are manually set in order to span a wide spectrum of time-varying local grouping patterns (see Figure below); more specifically:

- **bases 1, 2 and 5** exhibit a single, structural change of the partition between `t = 5` and `t = 6`;
- **bases 3 and 4** produce stable clusters over time (all countries separated for basis 3, and a single country isolated from the others for basis 4)
- **basis 6** is the most challenging regime, with all the countries changing group membership frequently across the whole time window.

``` r
n  <- 5 # number of countries
TT <- 10 # number of time points
p  <- 6 # number of spline bases

source("simulate_data.R") # load the main script
sim      <- simulate_scenario_4_1(n = n, TT = TT, p = p, seed = 290497)
clusters <- sim$clusters
```

The object `sim` contains the simulated log-mortality array `log_m` (countries × ages × periods), the true country-specific coefficients `beta` (countries × bases × periods), the cluster-specific values `theta`, the true memberships `clusters` and the B-spline basis matrix `S`.

The resulting cluster configuration is displayed below, reproducing Figure 3 in the manuscript.

``` r
## Common palette
gb <- paste0("#", c("fb4934", "b8bb26", "fabd2f", "83a598", "d3869b"))

df_pl <- reshape2::melt(clusters)
df_pl$value <- factor(df_pl$value)
df_pl$Var1  <- ordered(factor(gsub("Unit", "", df_pl$Var1)), levels = n:1)
df_pl$Var2  <- factor(df_pl$Var2)

clust_pl <- ggplot(df_pl) +
  geom_tile(aes(Var2, Var1, fill = value, width = 1, height = 0.95),
            linewidth = .5, col = "black", show.legend = FALSE, alpha = .8) +
  geom_text(aes(Var2, Var1, label = value), col = "white", size = 5,
            show.legend = FALSE) +
  facet_wrap(~L1) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey70")) +
  scale_x_discrete(expand = expansion(add = c(0.6))) +
  scale_y_discrete(expand = expansion(add = c(0.55))) +
  scale_fill_manual(values = gb) +
  xlab("Time") + ylab("Countries")

# plot(clust_pl)

ggsave(clust_pl, file = "img/sim_setting.png", width = 15, height = 7, dpi = 150)
```

![True cluster assignments in the simulation study](img/sim_setting.png)

## 2. Posterior computation

Posterior inference proceeds under the model proposed in Section 2 of the manuscript with diffuse hyperparameters `a = b = a_tau = b_tau = a_lambda = b_lambda = 1e-3`, `a_M = 2e-3`, `b_M = 1e-3` and `a_alpha = b_alpha = 1`. The entries of the GP covariance matrix on the coefficients’ means (eq (4) of the manuscript) are defined through a squared-exponential kernel with length scale `1.5`, while the mean vectors $\mu_j$ are elicited in a data-driven manner: for each period, the spline coefficients are first estimated via OLS under, and then smoothed over time via LOESS; refer to Section 4.1 of the manuscript for further details on this strategy.

The Gibbs sampler is run for `20000` iterations, discarding the first `10000` as a burn-in. The retained samples are then summarised into two quantities of interest:

- **the cluster memberships.** For each basis `j = 1, ..., 6` and period `t = 1, ..., 10`, a single point estimate $\hat{c}_{jt}$ is obtained by minimising the posterior expected variation of information (Wade & Ghahramani, 2018), using the function `minVI()` in the `mcclust.ext` package. The resulting partitions are stored in `coclust.RDS`, together with the true cluster allocations. In order to provide a comprehensive assessment for the clustering accuracy, the posterior mean of the percentage of pairs of countries that are correctly co-clustered is also computed.
- **the spline coefficients.** Point and $95\%$ interval estimates for the time- and country-specific coefficients are obtained as posterior means and quantile-base intervals, and stored in `beta_countries.RDS`

The chunk below runs the Gibbs sampler and saves the full output in the file `res_simstudy_scenario1.RDS` (currently gitignored). The outputs of interest are post-processed into the files `coclust.RDS` and `beta_countries.RDS`, as described below.

On a `M5` Mac-Air with an optimized `openblas` library (0.3.33) the running time is about 3 minutes.

``` r
out_MCMC <- run_model(Y = sim$log_m, ages = sim$ages,
  n_iter = 2e4, print_step = 500,
  name_save = "res_simstudy_scenario1.RDS")

burnin = 1:1e4

# Inference on partition (single point estimates, not used directly in the tutorial)
psm               <- compute_prob_coclust(import = FALSE, save = FALSE, draws = out_MCMC$res$labels, burnin = burnin)
ppe               <- compute_partition_point_est(psm = psm, draws = out_MCMC$res$labels, save = FALSE, burnin = burnin)

# Co-clustering summary
psm_acc <- summarise_clusters(draws = out_MCMC$res$labels, clusters_true = clusters, burnin = burnin)
coclust <- list(psm = psm, ppe = ppe, true = clusters, psm_acc = psm_acc)

# Inference on individual coefficients
beta_est       <- compute_beta_units(draws_beta = out_MCMC$res$beta, draws_labs = out_MCMC$res$labels, burnin = burnin)
beta_countries <- summarise_beta(beta_est, sim$beta)

# Save both objects
saveRDS(coclust, "coclust.RDS")
saveRDS(beta_countries, "beta_countries.RDS")
```

## 3. Results and performance

The co-clustering accuracies are reported below, as in Table 1 of the manuscript.

``` r
coclust <- readRDS("coclust.RDS")

knitr::kable(coclust$psm_acc, digits = 3, caption = "Posterior means of co-clustering accuracies for each combination (j, t).")
```

|         |     1 |     2 |     3 |     4 |     5 |     6 |     7 |     8 |     9 |    10 |
|:--------|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| Spline1 | 0.994 | 1.000 | 1.000 | 0.987 | 0.952 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| Spline2 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.991 | 1.000 | 1.000 | 1.000 | 0.998 |
| Spline3 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Spline4 | 0.999 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 |
| Spline5 | 0.999 | 0.999 | 1.000 | 0.998 | 0.992 | 0.998 | 1.000 | 1.000 | 1.000 | 0.997 |
| Spline6 | 0.942 | 0.917 | 0.965 | 0.964 | 0.984 | 0.981 | 0.979 | 0.976 | 0.988 | 0.990 |

Posterior means of co-clustering accuracies for each combination (j, t).

The figure below (Figure 4 in the paper) compares the estimated dynamic spline coefficients (lines) and their $95%$ credible intervals (ribbons) with the true ones (triangles), with one panel per spline basis and one colour per country.

``` r
beta_df <- readRDS("beta_countries.RDS")
beta_df$Unit   <- gsub("Unit", "Country ", beta_df$Unit)
beta_df$Spline <- gsub("Spline", "Spline ", beta_df$Spline)

beta_pl <- ggplot(beta_df) +
  geom_ribbon(aes(Year, ymin = Post_q025, ymax = Post_q975, fill = Unit),
              alpha = .25, show.legend = FALSE) +
  geom_line(aes(Year, Post_Mean, col = Unit), alpha = .9,
            show.legend = FALSE, linewidth = .9) +
  geom_point(aes(Year, Beta_true, col = Unit), size = 2, shape = "triangle") +
  facet_wrap(~Spline, scales = "free") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = gb) +
  scale_fill_manual(values = gb) +
  guides(color = guide_legend(title = NULL,
                              override.aes = list(shape = 22, size = 5, fill = gb))) +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey70"),
        legend.background = element_rect(fill = "grey", color = "black"),
        legend.key        = element_rect(fill = "white"),
        legend.position   = "bottom",
        legend.direction  = "horizontal") +
  xlab("Time") + ylab("Coefficients")
# plot(beta_pl)
ggsave(beta_pl, file = "img/sim_beta.png", width = 14, height = 7, dpi = 150)
```

![Estimated and true B-spline coefficients](img/sim_beta.png)
