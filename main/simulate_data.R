##########################################################################
## Data-generating mechanisms for the simulation studies in Section 4    ##
## of "Bayesian local clustering of age-period mortality surfaces        ##
## across multiple countries" (Romano, Aliverti and Durante)             ##
##########################################################################

## ----------------------------------------------------------------------
## B-spline bases
## ----------------------------------------------------------------------
## Quadratic B-spline bases g_1(x), ..., g_p(x) evaluated on the age grid.
## The bases are equally spaced over the range of `ages` and are the ones
## entering the expansion f_it(x) = sum_j beta_ijt g_j(x) in equation (2).

spline_basis <- function(ages, p = 6, degree = 2) {
  B <- splines::bs(ages, df = p, degree = degree, intercept = TRUE)
  B <- matrix(as.numeric(B), nrow = length(ages), ncol = p)
  rownames(B) <- ages
  colnames(B) <- paste0("Spline", seq_len(p))
  B
}

## ----------------------------------------------------------------------
## True local clustering structure of the first scenario (Section 4.1)
## ----------------------------------------------------------------------
## Cluster memberships c_jt are NOT drawn from the tRPM prior: they are set
## manually to span a wide spectrum of time-varying local grouping patterns.
## The output is a list of p = 6 matrices of dimension n x T (5 countries by
## 10 periods); entry [i, t] is the cluster label of country i at period t
## for the j-th spline basis.

cluster_setting_4_1 <- function() {

  ## bases 1-3
  b0 <- rbind(
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 2, 2, 2, 2, 2,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  2, 2, 2, 2, 2, 5, 5, 5, 5, 5,  2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    c(1, 1, 1, 1, 1, 3, 3, 3, 3, 3,  4, 4, 4, 4, 4, 2, 2, 2, 2, 2,  3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    c(2, 2, 2, 2, 2, 4, 4, 4, 4, 4,  3, 3, 3, 3, 3, 5, 5, 5, 5, 5,  4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
    c(2, 2, 2, 2, 2, 5, 5, 5, 5, 5,  4, 4, 4, 4, 4, 5, 5, 5, 5, 5,  5, 5, 5, 5, 5, 5, 5, 5, 5, 5))

  ## bases 4-6
  b1 <- rbind(
    c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  1, 1, 1, 1, 1, 3, 3, 3, 3, 3,  3, 3, 4, 4, 5, 5, 3, 3, 2, 2),
    c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  1, 1, 1, 1, 1, 4, 4, 4, 4, 4,  4, 4, 1, 1, 3, 3, 1, 1, 4, 4),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  2, 2, 2, 2, 2, 5, 5, 5, 5, 5,  3, 3, 5, 5, 5, 5, 5, 5, 2, 2),
    c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  1, 1, 1, 1, 1, 5, 5, 5, 5, 5,  4, 4, 4, 4, 3, 3, 5, 5, 1, 1),
    c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3, 3, 3, 3, 3, 5, 5, 5, 5, 5,  2, 2, 4, 4, 2, 2, 3, 3, 5, 5))

  list("Spline 1" = b0[, 1:10],  "Spline 2" = b0[, 11:20], "Spline 3" = b0[, 21:30],
       "Spline 4" = b1[, 1:10],  "Spline 5" = b1[, 11:20], "Spline 6" = b1[, 21:30])
}

## ----------------------------------------------------------------------
## Synthetic log-mortality rates for the first scenario (Section 4.1)
## ----------------------------------------------------------------------
## Data are generated exactly under model (1)-(4):
##
##   log m_ixt = f_it(x) + eps_ixt,           eps_ixt ~ N(0, sigma^2)
##   f_it(x)   = sum_{j=1}^p beta_ijt g_j(x)
##   beta_ijt  = theta_{c_ijt, j, t}
##   theta_kjt ~ N(gamma_jt, tau_j^2),        gamma_jt = gamma_j0 + slope * (t - 1)
##
## Arguments
##   clusters   list of p matrices (n x T) with the true memberships c_jt
##   ages       age grid X
##   sigma      country-specific residual sd (sigma_i = sigma for all i)
##   tau        vector of length p with the cluster-level sds tau_j
##   intercepts vector of length p with the intercepts gamma_j0 of the
##              parallel decreasing lines gamma_j
##   slope      common slope of the lines gamma_j
##
## Value: a list with the simulated log-mortality array `log_m` (n x |X| x T),
## the true country-specific coefficients `beta` (n x p x T), the cluster-
## specific values `theta`, the memberships `clusters` and the basis `B`.

simulate_scenario_4_1 <- function(clusters   = cluster_setting_4_1(),
                                  ages       = 0:100,
                                  sigma      = 0.05,
                                  tau        = c(0.05, 0.05, 0.05, 0.05, 0.10, 0.05),
                                  intercepts = c(-3.60, -5.82, -4.03, -5.05, -1.80, -0.73),
                                  slope      = -0.02) {

  p <- length(clusters)
  n <- nrow(clusters[[1]])
  TT <- ncol(clusters[[1]])
  stopifnot(length(tau) == p, length(intercepts) == p)

  B <- spline_basis(ages, p = p)

  ## higher-level mean trajectories gamma_j: parallel decreasing lines
  gamma <- outer(intercepts, slope * (seq_len(TT) - 1), "+")   # p x T

  ## cluster-specific coefficients theta_kjt and country-specific betas
  theta <- vector("list", p)
  beta  <- array(NA_real_, dim = c(n, p, TT),
                 dimnames = list(paste0("Unit", seq_len(n)),
                                 paste0("Spline", seq_len(p)),
                                 seq_len(TT)))

  for (j in seq_len(p)) {
    theta[[j]] <- vector("list", TT)
    for (t in seq_len(TT)) {
      lab <- clusters[[j]][, t]
      K   <- max(lab)
      th  <- rnorm(K, mean = gamma[j, t], sd = tau[j])
      theta[[j]][[t]] <- th
      beta[, j, t] <- th[lab]
    }
  }

  ## log-mortality surfaces plus Gaussian noise
  log_m <- array(NA_real_, dim = c(n, length(ages), TT),
                 dimnames = list(paste0("Unit", seq_len(n)), ages, seq_len(TT)))
  for (t in seq_len(TT)) {
    f <- beta[, , t] %*% t(B)                       # n x |X|
    log_m[, , t] <- f + rnorm(n * length(ages), 0, sigma)
  }

  list(log_m = log_m, beta = beta, theta = theta,
       clusters = clusters, gamma = gamma, B = B, ages = ages)
}
