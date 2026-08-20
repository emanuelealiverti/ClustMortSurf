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
##
## Arguments
##
##    ages: a vector of age values
##    p: the number of B-spline bases
##    knots: a vector of knots for the B-spline bases; if NULL, p equally spaced.
##    degree: the degree of the B-spline bases
##
## Value: a matrix of dimension length(ages) x p containing the B-spline bases 
## evaluated at the ages.

spline_basis <- function(ages, knots = NULL, p = NULL, degree = 2) {
  
  if (is.null(knots) & is.null(p)) {
    stop("Either knots or p must be specified.")
  }
  
  if (!is.null(knots) & !is.null(p)) {
    stop("Either knots or p must be specified, not both.")
  }
  
  require(splines2)
  
  S.temp <- bSpline(ages, knots = knots, degree = degree, intercept = TRUE)
  K <- ncol(S.temp)
  for (k in 1:K){S.temp[,k] <- S.temp[,k]/max(S.temp[,k])}
  
  # Add Dirac delta in 0 as first spline
  S <- rbind(c(1, rep(0, k)), cbind(rep(0, length(ages)), S.temp))
  
  return(S)
}



## ----------------------------------------------------------------------
## True local clustering structure of the first scenario (Section 4.1)
## ----------------------------------------------------------------------
## Cluster memberships c_jt are NOT drawn from the tRPM prior: they are set
## manually to span a wide spectrum of time-varying local grouping patterns.
## The output is a list of p = 6 matrices of dimension n x T (5 countries by
## 10 periods); entry [i, t] is the cluster label of country i at period t
## for the j-th spline basis.
##
## Value: a list of p matrices (n x T) with the true memberships c_jt.

cluster_setting_4_1 <- function(n, TT) {
  
  
  
  #1
  b1 = matrix(c(rep(1, TT/2), rep(1, TT/2),
                rep(1, TT/2), rep(1, TT/2), 
                rep(1, TT/2), rep(3, TT/2),
                rep(2, TT/2), rep(4, TT/2),
                rep(2, TT/2), rep(5, TT/2)),
              nrow = n, ncol = TT, byrow = T)
  #2
  b2 = matrix(c(rep(1, TT/2), rep(2, TT/2),
                rep(2, TT/2), rep(5, TT/2),
                rep(4, TT/2), rep(2, TT/2),
                rep(3, TT/2), rep(5, TT/2),
                rep(4, TT/2), rep(5, TT/2)),
              nrow = n, ncol = TT, byrow = T)
  #3
  b3 = matrix(c(rep(1, TT/2), rep(1, TT/2),
                rep(2, TT/2), rep(2, TT/2),
                rep(3, TT/2), rep(3, TT/2),
                rep(4, TT/2), rep(4, TT/2),
                rep(5, TT/2), rep(5, TT/2)),
              nrow = n, ncol = TT, byrow = T)
  #4
  b4 = matrix(c(rep(3, TT/2), rep(3, TT/2),
                rep(3, TT/2), rep(3, TT/2),
                rep(1, TT/2), rep(1, TT/2),
                rep(3, TT/2), rep(3, TT/2),
                rep(3, TT/2), rep(3, TT/2)),
              nrow = n, ncol = TT, byrow = T)
  #5
  b5 = matrix(c(rep(1, TT/2), rep(3, TT/2),
                rep(1, TT/2), rep(4, TT/2),
                rep(2, TT/2), rep(5, TT/2),
                rep(1, TT/2), rep(5, TT/2),
                rep(3, TT/2), rep(5, TT/2)),
              nrow = n, ncol = TT, byrow = T)
  #6
  b6 = matrix(c(rep(sample(1:n, n, TRUE), each = 2),
                rep(sample(1:n, n, TRUE), each = 2),
                rep(sample(1:n, n, TRUE), each = 2),
                rep(sample(1:n, n, TRUE), each = 2),
                rep(sample(1:n, n, TRUE), each = 2)),
              nrow = n, ncol = TT, byrow = T)
  
  dimnames(b1) <- dimnames(b2) <- dimnames(b3) <- 
    dimnames(b4) <- dimnames(b5) <- dimnames(b6) <- list(paste0("Unit", 1:5), 1:10)
  
  list("Spline 1" = b1,  "Spline 2" = b2, "Spline 3" = b3,
       "Spline 4" = b4,  "Spline 5" = b5, "Spline 6" = b6)
}



## ----------------------------------------------------------------------
## Synthetic log-mortality rates for the first scenario (Section 4.1)
## ----------------------------------------------------------------------
## Data are generated under:
##
##    log m_ixt = f_it(x) + eps_ixt,           eps_ixt ~ N(0, sigma^2)
##    f_it(x)   = sum_{j=1}^p beta_ijt g_j(x)
##    beta_ijt  = theta_{c_ijt, j, t}
##    theta_kjt ~ N(phi_kjt, tau_j^2),       phi_kjt = int_{kjt} + slope * (t - 1)
##
## Arguments
##    clusters        list of p matrices (n x T) with the true memberships c_jt
##    ages            age grid X
##    sigma           country-specific residual sd (sigma_i = sigma for all i)
##    delta           vector of length p with the cluster-level sds delta_j
##    min_intercepts  vector of length p with the intercepts gamma_j0 of the
##                    parallel decreasing lines gamma_j
##    max_intercepts  vector of length p with the intercepts gamma_j0 of the
##                    parallel decreasing lines mu_j
##    slope           common slope of the lines mu_j
##
## Value: a list with the simulated log-mortality array `log_m` (n x |X| x T),
## the true country-specific coefficients `beta` (n x p x T), the cluster-
## specific values `theta`, the memberships `clusters` and the basis `S`.

simulate_scenario_4_1 <- function(n               = 5, 
                                  TT              = 10, 
                                  p               = 6, 
                                  ages            = 0:100,
                                  knots           = c(20, 40),
                                  min_intercepts  = c(-4.4, -7.1, -5, -5.4, -2.5, -1.5),
                                  max_intercepts  = c(-2, -4.5, -3, -3.8, -1, 0),
                                  sigma           = c(0.05, 0.05, 0.05, 0.05, 0.05),
                                  delta           = c(0.05, 0.10, 0.05, 0.05, 0.05, 0.05),
                                  slope           = -0.02,
                                  seed            = 290497) {
  
  require(mvtnorm)
  
  set.seed(seed)
  
  stopifnot(length(delta) == p, length(min_intercepts) == p, length(max_intercepts) == p)
  
  if (!is.null(knots)) {
    S <- spline_basis(ages[-1], knots = knots)
  } else {
    S <- spline_basis(ages[-1], p = p)
  }
  
  # Correlation matrix of GP
  grid_years <- expand.grid(1:TT, 1:TT)
  val_ker <- apply(grid_years, 1, function(x) sq_exp_ker(x[1] - x[2], 3/2))
  Sigma <- matrix(val_ker, TT, TT, byrow = F)
  
  
  ## higher-level mean trajectories mu_j: parallel decreasing lines
  intercepts <- matrix(NA, nrow = n, ncol = p)
  for (j in seq_len(p)){
    intercepts[ , j] = seq(min_intercepts[j], max_intercepts[j], length.out = n)
  }
  phi <- array(NA_real_, dim = c(n, p, TT),
               dimnames = list(paste0("Cluster", seq_len(n)),
                               paste0("Spline", seq_len(p)),
                               seq_len(TT)))
  for (j in seq_len(p)) {
    for (k in seq_len(n)){
      phi[k, j, ] <- intercepts[k, j] + slope * seq_len(TT)
    }
  }
  
  
  ## cluster-specific coefficients beta_kjt and country-specific betas
  theta <- array(NA_real_, dim = c(n, p, TT),
                 dimnames = list(paste0("Cluster", seq_len(n)),
                                 paste0("Spline", seq_len(p)),
                                 seq_len(TT)))
  beta  <- array(NA_real_, dim = c(n, p, TT),
                 dimnames = list(paste0("Unit", seq_len(n)),
                                 paste0("Spline", seq_len(p)),
                                 seq_len(TT)))
  
  for (j in seq_len(p)) {
    for (k in seq_len(n)){
      th <- drop(rmvnorm(1, phi[k, j, ], delta[j]^2*Sigma))
      theta[k, j, ] <- th
    }
  }
  
  clusters <- cluster_setting_4_1(n = n, TT = TT)
  for (j in seq_len(p)){
    for (t in seq_len(TT)) {
      lab <- clusters[[j]][, t]
      beta[, j, t] <- theta[lab, j, t]
    }
  }
  
  ## log-mortality surfaces plus Gaussian noise
  log_m <- array(NA, dim = c(n, length(ages), TT),
                 dimnames = list(paste0("Unit", seq_len(n)), ages, seq_len(TT)))
  # for (t in seq_len(TT)) {
  #   f <- beta[, , t] %*% t(S)                       # n x |X|
  #   log_m[, , t] <- f + rnorm(n * length(ages), 0, sigma[i])
  # }
  for (i in seq_len(n)) {
    f <- t(beta[i, , ]) %*% t(S)                       # n x |X|
    log_m[i, , ] <- t(f + rnorm(length(ages)*TT, 0, sigma[i]))
  }
  
  list(log_m = log_m, beta = beta, theta = theta,
       clusters = clusters, S = S, ages = ages, seed = seed)
}
