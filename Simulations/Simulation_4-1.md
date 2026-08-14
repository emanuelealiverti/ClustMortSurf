# Simulation study — first scenario (Section 4.1)

This tutorial reproduces the **first simulation scenario** of the paper
[*Bayesian local clustering of age-period mortality surfaces across multiple countries*](https://arxiv.org/abs/2504.05240),
in which synthetic log-mortality rates are generated **under the proposed model**.
 The goal of the study is to assess to what extent the proposed formulation is able to recover (i) complex local grouping structures among countries that vary across ages and periods, and (ii) the associated cluster-specific coefficients.

The analysis proceeds in three steps:

1. [simulating the data](#1-data-generating-mechanism) and displaying the [true local clustering structure](#2-true-local-clustering-structure);
2. [running the Gibbs sampler](#3-posterior-computation) and storing its two posterior summaries of interest, i.e. the estimated cluster memberships in `coclust_df.RDS` and the posterior means of the spline coefficients in `beta_df.RDS`;
3. [assessing the recovery of the local grouping structures](#4-co-clustering);
4. [comparing estimated and true trajectories](#5-estimated-versus-true-coefficients) of the country-specific coefficients.

The code is meant to be run from within this folder, and the core routines are sourced from [`main/`](../main). 

```r
library(ggplot2)

source("../main/simulate_data.R")

set.seed(4321)

## Common palette for clusters and spline bases
gb <- paste0("#", c("fb4934", "b8bb26", "fabd2f", "83a598", "d3869b", "8ec07c"))
```

---

## 1. Data-generating mechanism

We simulate log-mortality rates for `n = 5` synthetic countries, observed over `T = 10` periods and
for the ages `X = {0, 1, ..., 100}`. Data are generated under equations (1)–(3) of the paper and splines
coefficient distribution centered around parallel lines:

```
log m_ixt = f_it(x) + eps_ixt,          eps_ixt ~ N(0, sigma_i^2)
f_it(x)   = sum_{j=1}^p beta_ijt g_j(x)
beta_ijt  = theta_{c_ijt, j, t}
theta_kjt ~ N(phi_kjt, delta_j^2)        phi_kjt = intercept_{kjt} + slope * (t - 1)
```

so that the age structure of mortality is described by a B-spline expansion whose country-specific
coefficients `beta_ijt` are **tied** whenever two countries belong to the same cluster at the pair
(basis `j`, period `t`). 

Note that `(phi_kj1, ..., phi_kjT)` for different clusters `k` are parallel lines in time:
the paper fixes their common slope (`-0.02`), while their intercepts are chosen so that 
the resulting surfaces reproduce realistic levels of log-mortality across the age range. 
The default intercepts in `simulate_scenario_4_1()` are the ones used to produce the data analysed below.

```r
clusters <- cluster_setting_4_1()
sim      <- simulate_scenario_4_1(clusters = clusters, ages = 0:100)

str(sim$log_m)
#  num [1:5, 1:101, 1:10] -3.6 -3.56 -3.71 -3.6 -3.62 ...
#  - attr(*, "dimnames")=List of 3
#   ..$ : chr [1:5] "Unit1" "Unit2" "Unit3" "Unit4" ...
#   ..$ : chr [1:101] "0" "1" "2" "3" ...
#   ..$ : chr [1:10] "1" "2" "3" "4" ...
```

The object `sim` contains the simulated log-mortality array `log_m` (countries × ages × periods), the true country-specific coefficients `beta` (countries × bases × periods), the cluster-specific values `theta`, the true memberships `clusters` and the B-spline basis matrix `B`.

## 2. True local clustering structure

The cluster memberships `c_jt` are manually set as in `cluster_setting_4_1()` — in order to span a wide spectrum of time-varying local grouping patterns, and thereby stress-test the model in regimes that are progressively harder to recover:

- **bases 3 and 4** produce stable clusters over time (all countries separated for basis 3, and a single country isolated from the others for basis 4), albeit with rather different co-clustering patterns;
- **bases 1, 2 and 5** exhibit a single, structural change of the partition halfway through the observational window (between `t = 5` and `t = 6`);
- **basis 6** is the most challenging regime, with all the countries changing group membership frequently across the whole time window.

The resulting configuration is displayed below, (see also Figure 3 in the paper).

```r
df_pl <- reshape2::melt(clusters)
df_pl$value <- factor(df_pl$value)
df_pl$Var1  <- ordered(factor(df_pl$Var1), levels = 5:1)
df_pl$Var2  <- factor(df_pl$Var2)

plf <- ggplot(df_pl) +
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

ggsave(plf, file = "img/sim_setting.png", width = 15, height = 7, dpi = 150)
```

![True cluster assignments in the simulation study](img/sim_setting.png)

## 3. Posterior computation

Posterior inference proceeds under the model proposed in Section 2 of the manuscript with diffuse hyperparameters
`a = b = a_tau = b_tau = a_lambda = b_lambda = 1e-3`, `a_M = 2e-3`, `b_M = 1e-3` and
`a_alpha = b_alpha = 1`. The entries of the GP covariance matrix `Lambda` are defined through a
squared-exponential kernel with length scale 1.5, i.e. `Lambda[t,t'] = exp(-0.5 (t - t')^2 / 1.5^2)`,
while the mean vectors `mu_j` are elicited in a data-driven manner: for each period, the spline
coefficients are first estimated via OLS under (1)–(2), and then smoothed over time via LOESS.

The Gibbs sampler is run for 20,000 iterations, discarding the first 10,000 as a conservative
burn-in. Traceplots and autocorrelation plots indicate satisfactory mixing of the chains. The 10,000
retained samples are then summarised, as detailed in Section 3.2 of the manuscript, into the two
quantities of interest.

- **the cluster memberships.** For each basis `j = 1, ..., 6` and period `t = 1, ..., 10`, the
  samples of `c_jt` are first summarised into the posterior co-clustering matrix `P_jt`, whose
  generic entry is the proportion of samples in which countries `i` and `i'` share the same group.
  A single point estimate `c_jt` is then obtained by minimising the posterior expected variation of
  information (Wade & Ghahramani, 2018) — `minVI()` in the `mcclust.ext` package — using `P_jt` as
  input. The result is stored in
  `coclust_df.RDS`.
- **the spline coefficients.** Point estimates for the country-specific trajectories
  `beta_ijt = theta_{c_ijt, j, t}` are obtained as posterior means and stored in `beta_df.RDS` alongside the values used to generate the data.

```r
## source("../main/gibbs_tRPM.R")

## fit         <- gibbs_tRPM(log_m = sim$log_m, ages = sim$ages, p = 6,
##                           n_iter = 20000, burn_in = 10000, ...)

## coclust_df  <- summarise_clusters(fit, clusters_true = clusters)   ## via minVI
## beta_df     <- summarise_beta(fit, beta_true = sim$beta)

## saveRDS(coclust_df, "coclust_df.RDS")
## saveRDS(beta_df,    "beta_df.RDS")
```

> **Note.** To be added

## 4. Co-clustering

*To be completed once `coclust_df.RDS` is available.* This section will compare the estimated
partitions with the true ones displayed in Section 2, and report the co-clustering accuracies of
Table 1 of the paper, i.e. the posterior mean, for each pair `(j, t)`, of the percentage of pairs of
countries that are correctly co-clustered.

## 5. Estimated versus true coefficients

The file `beta_df.RDS` stores, for every country (`Unit`), spline basis (`Spline`) and period
(`Year`), the posterior mean of the country-specific coefficient (`Beta_est`) alongside the value
used to generate the data (`Beta_true`).

```r
beta_df <- readRDS("beta_df.RDS")
head(beta_df)
#   Spline  Unit  Year Beta_est Beta_true
#   Spline1 Unit1     1    -4.38     -4.43
#   Spline1 Unit1     2    -4.41     -4.41
#   Spline1 Unit1     3    -4.41     -4.45
#   Spline1 Unit1     4    -4.51     -4.52
#   Spline1 Unit1     5    -4.54     -4.55
#   Spline1 Unit1     6    -4.53     -4.55
```

The figure below (Figure 4 in the paper) compares the estimated trajectories (lines) with the true
ones (triangles), with one panel per country and one colour per spline basis.

```r
beta_df$Unit   <- gsub("Unit", "Country ", beta_df$Unit)
beta_df$Spline <- gsub("Spline", "Spline ", beta_df$Spline)

beta_pl <- ggplot(beta_df) +
  geom_line(aes(Year, Beta_est, col = Spline), alpha = .9,
            show.legend = FALSE, linewidth = .9) +
  geom_point(aes(Year, Beta_true, col = Spline), size = 2, shape = "triangle") +
  facet_wrap(~Unit, scales = "free") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = gb) +
  guides(color = guide_legend(title = NULL,
                              override.aes = list(shape = 22, size = 5, fill = gb))) +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey70"),
        legend.background = element_rect(fill = "grey", color = "black"),
        legend.key        = element_rect(fill = "white"),
        legend.position   = "right",
        legend.direction  = "horizontal") +
  xlab("Time") + ylab("Coefficients")

## the legend is moved into the empty panel of the facet grid
beta_pl <- lemon::reposition_legend(beta_pl, "center", panel = "panel-3-2",
                                    plot = FALSE)

ggsave(beta_pl, file = "img/sim_beta.png", width = 14, height = 7, dpi = 150)
```

![Estimated and true B-spline coefficients](img/sim_beta.png)

