###################################################################################
## Implementation of the full-conditional updates and the Gibbs sampler          ##
## routine for "Bayesian local clustering of age-period mortality                ##
## surface across multiple countries" (Romanò, Aliverti and Durante, 2025+)      ##
###################################################################################



## ----------------------------------------------------------------------
## Update temporal transition indicator (gamma)
## ----------------------------------------------------------------------
## Updates gamma_{it}, which indicates whether unit i stays in the same cluster
## between time t-1 and time t.
##
## Arguments
##    i           index of the observation unit to update
##    gamma       vector of gamma_{t} transition indicators
##    alpha_t     tRPM probability parameter at time t
##    lab_t       cluster labels vector at time t
##    lab_tm1     cluster labels vector at time t-1
##    M           CRP concentration parameter
##
## Value: binary integer (0 or 1) representing the updated gamma_{it}.

up_gamma_i <- function(i, gamma, alpha_t, lab_t, lab_tm1, M){
  
  # R_t: set of units fixed from t-1 to t, including unit i
  # R_tmi: R_t \ {i}
  # R_tpi: R_tmi U {i}
  R_t <- which(gamma == 1)
  R_tmi <- R_t[R_t != i]
  R_tpi <- c(R_tmi, i)
  
  # Reduced partition labels
  lab_t.R <- rep(-1, length(lab_t))
  lab_t.R[R_tpi] <- lab_t[R_tpi]
  
  # Check if co-clustering matches between t-1 and t on the reduced partition
  adj_tm1 <- outer(lab_tm1[R_tpi], lab_tm1[R_tpi], function(x, y) as.integer(x == y))
  adj_t <- outer(lab_t[R_tpi], lab_t[R_tpi], function(x, y) as.integer(x == y))
  check <- all(adj_tm1 == adj_t)
  
  # If check is FALSE, computation is omitted as probability will be 0
  if (check){
    # Find partition set containing unit i
    j <- lab_t.R[i]
    n.R <- sum(lab_t.R > 0) - 1
    n_j <- sum(lab_t.R == j)
    
    if (n.R == 0) {
      ratio <- 1
    } else {
      if (n_j == 1) {
        # Unit i forms a new singleton cluster
        ratio <- M / (n.R + M)
      } else {
        # Unit i belongs to an existing cluster
        ratio <- (n_j - 1) / (n.R + M)
      }
    }
    
    prob <- alpha_t / (alpha_t + (1 - alpha_t) * ratio) * check
    gamma_it <- as.integer(runif(1) < prob)
    
  } else {
    gamma_it <- 0
  }
  
  return(gamma_it)
}



## ----------------------------------------------------------------------
## Update cluster membership label
## ----------------------------------------------------------------------
## Updates the cluster assignment for a specific unit, spline coefficient,
## and time point.
##
## Arguments
##    i               index of the observation unit
##    j               index of the spline coefficient
##    Y_it            vector of responses Y_{ixt} across ages x
##    beta_i          vector of spline coefficients for unit i at time t
##    beta_cluster    vector of spline coefficients for current clusters at time t
##    sigma_i         standard deviation of observation i
##    lab_t           cluster labels vector at time t
##    lab_tp1         cluster labels vector at time t+1
##    gamma_tp1       gamma_{t+1} transition indicators (or 'last time' for final T)
##    newclustervalue sampled value for potential new cluster coefficient
##    spline_basis    matrix of evaluated spline basis functions
##    M               CRP concentration parameter
##
## Value: integer label indicating the updated cluster assignment c_{it}.

up_label_i <- function(i, j,
                       Y_it, beta_i, beta_cluster, sigma_i,
                       lab_t, lab_tp1, gamma_tp1,
                       newclustervalue, spline_basis, M) {
  
  # Sample size and cluster unit i belongs to
  n <- length(lab_t)
  k <- lab_t[i]
  
  # Partition excluding unit i
  lab_tmi <- lab_t
  lab_tmi[i] <- -1
  
  # Active clusters with at least one observation
  non_empty_clust <- whichclusters <- unique(lab_tmi[lab_tmi > 0])
  
  # Include candidate new cluster index
  whichclusters <- c(whichclusters, max(whichclusters) + 1L)
  
  # Partial Gaussian means without j-th spline
  means <- drop(spline_basis[, -j] %*% beta_i[-j])
  
  # Initialize storing objects
  beta_list <- c()
  prob_cit_list <- c()
  check_list <- c()
  
  spl_bas_j <- spline_basis[, j]
  
  is.last_time <- identical(gamma_tp1, 'last time')
  # If it is not last time compute co-clustering for compatibility
  if (!(is.last_time)){
    R_tp1 <- which(gamma_tp1 == 1)
    R_tp1_mi <- R_tp1[R_tp1 != i]
    adj_tp1 <- outer(lab_tp1[R_tp1], lab_tp1[R_tp1], function(x, y) as.integer(x == y))
  }
  
  for (h in whichclusters){
    # Verify compatibility of potential allocation of i to cluster h (for every 
    # potential allocation)
    if (is.last_time){
      check <- TRUE
    } else {
      lab_t.h <- lab_t
      lab_t.h[i] <- h
      adj_t.h <- outer(lab_t.h[R_tp1], lab_t.h[R_tp1], function(x, y) as.integer(x == y))
      check <- all(adj_t.h == adj_tp1)
    }
    
    # Compute prior contribution to probability of allocations.
    # If allocation to h is compatible use CRP predictive. Otherwise it is 0.
    if (check){
      if (h %in% non_empty_clust){
        beta <- beta_cluster[h]
        prob_cit <- sum(lab_tmi == h)
      } else {
        if (sum(lab_t == k) == 1){
          beta <- beta_cluster[k]
        } else{
          beta <- newclustervalue
        }
        prob_cit <- M
      }
    } else {
      prob_cit <- 0
      beta <- 0
    }
    
    beta_list <- c(beta_list, beta)
    prob_cit_list <- c(prob_cit_list, prob_cit)
    check_list <- c(check_list, check)
  }
  
  # Compute likelihood contribution to full-conditional
  logpnorm_list <- colSums(-0.5 * log(2 * pi * sigma_i^2) - 
                             0.5 / sigma_i^2 * (Y_it - means - outer(spl_bas_j, beta_list))^2)
  logp <- log(prob_cit_list) + logpnorm_list + log(check_list)
  
  # Sample through log-prob with max-Gumbel trick
  gumbel <- -log(-log(runif(length(logp))))
  lll <- logp + gumbel
  c_it <- whichclusters[which.max(lll)]
  
  # If the observation is assigned to a "new" cluster and it was in a singleton
  # then its new label will stay as it was before.
  if (c_it == max(whichclusters) && (sum(lab_t == k) == 1)){
    c_it = k
  }
  
  return(c_it)
}



## ----------------------------------------------------------------------
## Update spline cluster coefficients (beta)
## ----------------------------------------------------------------------
## Updates cluster-specific spline coefficient via Gaussian-Gaussian conjugacy.
##
## Arguments
##    j             index of the coefficient/spline to update
##    k             index of the cluster to update
##    t             time point index
##    Y_t           matrix with Y_{ixt} observations for time t (for every i, x)
##    beta_t        matrix of cluster spline coefficients at time t (for every spline)
##    phi_jt        prior mean for coefficient j at time t
##    delta_j       prior standard deviation for coefficient j
##    sigma_vec     vector of observation standard deviations
##    labels_t      matrix of cluster labels at time t
##    spline_basis  matrix of evaluated spline basis functions
##
## Value: scalar value containing the updated beta_{jkt}.

up_beta <- function(j, k, t,
                    Y_t, beta_t, phi_jt, delta_j, sigma_vec,
                    labels_t, spline_basis){
  
  # Compute sample size and numer of splines
  n <- length(sigma_vec)
  p <- ncol(spline_basis)
  
  # Get units i that are in cluster k for spline j
  obs <- which(labels_t[, j] == k)
  
  # Get spline coefficients (excluded j-th spline) for units in cluster k
  beta_minusj_units_clusterk <- vapply(1:p, 
                                       function(x) beta_t[labels_t[, x], x],
                                       FUN.VALUE = double(n))[obs, -j]
  
  # Compute partial residuals r_{itx(j)}
  r_t <- Y_t[obs, ] - beta_minusj_units_clusterk %*% t(spline_basis[, -j])
  r_t <- drop(r_t)
  
  sigma_clusterk <- sigma_vec[obs]
  
  # Conjugate Gaussian update computations
  var.post <- (1 / delta_j^2 + sum(1 / sigma_clusterk^2) * sum(spline_basis[, j]^2))^(-1)
  mean.post <- var.post * (sum(t(t(r_t) * spline_basis[, j]) / sigma_clusterk^2) + phi_jt / delta_j^2)
  
  beta_updated <- rnorm(1, mean = mean.post, sd = sqrt(var.post))
  
  return(beta_updated)
}



## ----------------------------------------------------------------------
## Update tRPM probability parameter (alpha)
## ----------------------------------------------------------------------
## Sample updated alpha parameter from conjugate Beta posterior.
##
## Arguments
##    sum_data      sum of gamma transition indicators at current time
##    n             total number of units
##    priorshape1   first parameter of Beta prior
##    priorshape2   second parameter of Beta prior
##
## Value: updated alpha parameter sampled from Beta posterior.

up_alpha_j <- function(sum_data,
                       n,
                       priorshape1,
                       priorshape2){
  
  out <- rbeta(1, 
               priorshape1 + sum_data, 
               priorshape2 + n - sum_data)
  
  return(out)
}



## ----------------------------------------------------------------------
## Update CRP concentration parameter (M_j)
## ----------------------------------------------------------------------
## Function to update the Dirichlet Process / CRP concentration parameter M_j
## for a specific spline component j.
##
## Arguments
##    gamma_j     matrix (n x T_final) of cluster assignment indicator variables 
##                for the j-th spline component
##    labels_j    matrix (n x T_final) of cluster labels for the j-th spline component
##    M_j         current scalar value of the CRP concentration parameter M_j
##    a_M         shape hyperparameter for Gamma prior on M_j
##    b_M         rate hyperparameter for Gamma prior on M_j
##    T_final     total number of time steps
##    n           total number of observation units (e.g., countries)
##
## Value: a single scalar representing the updated M_j draw.

up_M <- function(gamma_j, labels_j, M_j, 
                 a_M, b_M, T_final, n) {
  
  # Identify constrained units
  constr <- apply(gamma_j, 2, function(x) which(x == 1), simplify = F)
  n_constr <- vapply(constr, length, FUN.VALUE = integer(1))
  
  # Sample auxiliary beta variables eta_t
  eta_temp <- vapply(1:T_final, function(t) {
    if (n == n_constr[t]) {
      1.0
    } else {
      rbeta(1, M_j + n_constr[t], n - n_constr[t])
    }
  }, FUN.VALUE = double(1))
  
  # Total number of active clusters at each time step
  k_temp <- apply(labels_j, 2, function(x) length(unique(x)))
  
  # Number of distinct constrained clusters
  k_constr_temp <- numeric(T_final)
  k_constr_temp[1] <- 0
  for (t in 2:T_final) {
    if (n_constr[t] > 0) {
      k_constr_temp[t] <- length(unique(labels_j[constr[[t]], t]))
    } else {
      k_constr_temp[t] <- 0
    }
  }
  
  # Number of new clusters
  deltak_temp <- k_temp - k_constr_temp
  
  # Draw from Gamma distribution
  M_new <- rgamma(1, 
                  shape = a_M + sum(deltak_temp), 
                  rate  = b_M - sum(log(eta_temp)))
  
  return(M_new)
}



## ----------------------------------------------------------------------
## Run Gibbs Sampler
## ----------------------------------------------------------------------
## Main wrapper function to initialize and run the MCMC Gibbs sampler for the 
## temporal Random Partition Model (tRPM) with B-splines.
##
## Arguments
##    Y           list of observation matrices, one per country/unit, 
##                with dimensions (time_points x ages)
##    seed        random seed for reproducibility
##    a_sigma     shape hyperparameter for inverse-gamma prior on observation variance
##    b_sigma     rate hyperparameter for inverse-gamma prior on observation variance
##    a_alpha     shape1 hyperparameter for Beta prior on alpha
##    b_alpha     shape2 hyperparameter for Beta prior on alpha
##    a_delta     shape hyperparameter for inverse-gamma prior on delta_j^2
##    b_delta     rate hyperparameter for inverse-gamma prior on delta_j^2
##    a_omega     shape hyperparameter for inverse-gamma prior on omega_j^2
##    b_omega     rate hyperparameter for inverse-gamma prior on omega_j^2
##    a_M         shape hyperparameter for Gamma prior on CRP concentration parameter M
##    b_M         rate hyperparameter for Gamma prior on CRP concentration parameter M
##    l           smoothness parameter of the squared exponential kernel
##    n_iter      total number of MCMC iterations
##    print_step  step frequency for printing MCMC progress
##    path_save   destination directory path for output saving
##    name_save   filename for saving the output list as an .RDS object
##
## Value: a list containing posterior MCMC chains ('res'), prior hyperparameters
##        ('parameters'), input data ('data'), and fitted metadata ('other').

run_model <- function(Y,
                      ages,
                      seed = 4238,
                      a_sigma = 0.001, b_sigma = 0.001,
                      a_alpha = 1, b_alpha = 1,
                      a_delta = 0.001, b_delta = 0.001,
                      a_omega = 0.001, b_omega = 0.001,
                      a_M = 0.002, b_M = 0.001,
                      l = 3/2,
                      n_iter = 5000, print_step = 500,
                      path_save = "",
                      name_save = paste(paste("res", Sys.Date(), sep = "_"), ".RDS", sep = "")){
  
  # Transform input data into list of matrices (one per country/unit)
  Y <- apply(Y, 1, function(x) t(x), simplify = FALSE)
  
  library(splines2)
  set.seed(seed)
  
  # ----------------------------------------------------------------------
  # B-Spline basis matrix setup
  # ----------------------------------------------------------------------
  ages <- 0:100
  ages_no0 <- ages[ages!= 0]
  
  # Construct spline basis functions over ages
  S <- spline_basis(ages = ages_no0, knots = c(20, 40))
  
  # ----------------------------------------------------------------------
  # Model Dimensions and Covariance Matrix Initialization
  # ----------------------------------------------------------------------
  T_final <- nrow(Y[[1]]); cat("Number of considered years:", T_final, "\n")
  n <- length(Y); cat("Number of considered countries:", n, "\n")
  p <- ncol(S); cat("Number of considered splines:", p, "\n")
  
  # Compute Gaussian Process correlation matrix
  SIGMA <- outer(1:T_final, 1:T_final, function(x1, x2) sq_exp_ker(x1 - x2, l = l))
  cholSIG <- chol(SIGMA)
  SIGMAinv <- chol2inv(cholSIG)
  cholSIGinv <- solve(cholSIG)
  
  # ----------------------------------------------------------------------
  # Estimate Prior Mean-Level for Gaussian Processes via Regressions
  # ----------------------------------------------------------------------
  GP_means_tmp <- matrix(NA, nrow = p, ncol = T_final)
  cat("Start estimate of GP means \n")
  for (t in 1:T_final){
    cat(t, "")
    data_tmp <- sapply(Y, function(x) x[t, ])
    melt_tmp <- reshape2::melt(data_tmp)
    
    df <- data.frame(Y = melt_tmp$value, 
                     matrix(t(S), nrow = nrow(melt_tmp), ncol = p, byrow = TRUE))
    mod <- lm(Y ~ -1 + ., data = df)
    
    GP_means_tmp[, t] <- mod$coefficients
  }
  
  # Smooth the overall temporal trend using LOESS
  GP_means <- list()
  for (j in 1:p){
    df <- data.frame(x = 1:T_final, y = GP_means_tmp[j, ])
    GP_means[[j]] <- loess(y ~ x, data = df)$fitted
  }
  cat("\nEnd estimate of GP means \n\n")
  
  # ----------------------------------------------------------------------
  # Storage Objects for MCMC Output
  # ----------------------------------------------------------------------
  gamma_res <- labels_res <- replicate(p, array(NA, dim = c(n, T_final, n_iter)), simplify = FALSE)
  beta_res <- replicate(p, array(NA, dim = c(n, T_final, n_iter)), simplify = FALSE)
  phi_res <- replicate(p, array(NA, dim = c(T_final, n_iter)), simplify = FALSE)
  delta_res <- omega_res <- replicate(p, array(NA, dim = c(n_iter)), simplify = FALSE)
  alpha_res <- replicate(p, array(NA, dim = c(n_iter)), simplify = FALSE)
  sigma_res <- array(NA, dim = c(n, n_iter))
  M_res <- replicate(p, array(NA, dim = c(n_iter)), simplify = FALSE)
  
  # Temporary state structures
  gamma_temp <- labels_temp <- beta_temp <- replicate(p, array(NA, dim = c(n, T_final)), simplify = FALSE)
  phi_temp <- replicate(p, array(NA, dim = c(T_final)), simplify = FALSE)
  delta_temp <- omega_temp <- replicate(p, list(NA))
  alpha_temp <- replicate(p, list(NA))
  sigma_temp <- array(NA, dim = c(n))
  M_temp <- replicate(p, list(NA))
  
  # ----------------------------------------------------------------------
  # Initialize MCMC Chain
  # ----------------------------------------------------------------------
  sigma_res[, 1] <- sigma_temp <- sapply(Y, sd)
  
  cat("Start initialization of values for first iteration of Gibbs Sampler \n")
  for (j in 1:p){
    cat(j, "")
    
    omega_res[[j]][1] <- omega_temp[[j]] <- 1
    delta_res[[j]][1] <- delta_temp[[j]] <- 1
    
    phi_res[[j]][, 1] <- phi_temp[[j]] <- drop(rmvn(1,
                                                    mu = GP_means[[j]],
                                                    sigma = omega_temp[[j]] * cholSIG,
                                                    isChol = TRUE))
    
    alpha_res[[j]][1] <- alpha_temp[[j]] <- rbeta(1, a_alpha, b_alpha)
    gamma_res[[j]][,, 1] <- gamma_temp[[j]][,] <- 0
    M_temp[[j]] <- ceiling(n / 5)
    
    for (t in 1:T_final){
      lab <- rho2lab(rCRP(n, M_temp[[j]]))
      labels_temp[[j]][, t] <- labels_res[[j]][, t, 1] <- lab
      
      beta_res[[j]][1:max(lab), t, 1] <- beta_temp[[j]][1:max(lab), t] <- 
        rnorm(max(lab), mean = phi_temp[[j]][t], sd = delta_temp[[j]])
    }
  }
  cat("\nEnd initialization of values for first iteration of Gibbs Sampler \n\n")
  
  # ----------------------------------------------------------------------
  # Main Gibbs Loop
  # ----------------------------------------------------------------------
  cat("Start Gibbs Sampler -", format(Sys.time(), "%H:%M:%S"), "\n")
  inizio <- last_time <- Sys.time()
  
  for (d in 2:n_iter){
    
    for (j in 1:p){
      
      for (t in 1:T_final){
        
        # --- Update Gamma (transition indicators) ---
        for (i in 1:n){
          # The step for t = 1 can also be skipped because it is initialized as
          # a vector of 0 and it never changes.
          if (t == 1){
            gamma_temp[[j]][i, 1] <- 0
          } else {
            gamma_temp[[j]][i, t] <- up_gamma_i(i = i, 
                                                gamma = gamma_temp[[j]][, t], 
                                                alpha_t = alpha_temp[[j]],
                                                lab_t = labels_temp[[j]][, t],
                                                lab_tm1 = labels_temp[[j]][, t - 1],
                                                M = M_temp[[j]])
          }
        }
        
        # --- Update Cluster Labels ---
        # Get current iteration values of beta at time t for every country, spline
        beta_units <- vapply(1:p, 
                             function(x) beta_temp[[x]][labels_temp[[x]][, t], t],
                             FUN.VALUE = double(n))
        
        # Select units that can be re-assigned
        can_move <- which(gamma_temp[[j]][, t] == 0)
        
        for (i in can_move){
          
          # Sample value for potential new cluster
          newclustervalue <- rnorm(1, phi_temp[[j]][t], delta_temp[[j]])
          
          # If t = t+1 gamma_tp1 does not actually exist
          .gamma_tp1 <- if (t == T_final) 'last time' else gamma_temp[[j]][, t + 1]
          
          labels_temp[[j]][i, t] <- up_label_i(i = i, 
                                               j = j,
                                               Y_it = Y[[i]][t, ],
                                               beta_i = beta_units[i, ],
                                               beta_cluster = beta_temp[[j]][, t],
                                               sigma_i = sigma_temp[i],
                                               lab_t = labels_temp[[j]][, t],
                                               lab_tp1 = labels_temp[[j]][, t + 1], 
                                               gamma_tp1 = .gamma_tp1,
                                               newclustervalue = newclustervalue,
                                               spline_basis = S,
                                               M = M_temp[[j]])
          
          # Standardize labels sequentially without gaps
          no_gap <- dplyr::dense_rank(labels_temp[[j]][, t])
          labels_new <- order(unique(no_gap))[no_gap]
          
          # Reorder cluster spline coefficients to match label ordering
          if (all(labels_temp[[j]][i, t] > labels_temp[[j]][-i, t])){
            if (labels_temp[[j]][i, t] > n){
              beta_extended <- c(beta_temp[[j]][, t], newclustervalue)
              sort_beta_temp <- unique(labels_temp[[j]][, t])
              beta_new <- beta_extended[sort_beta_temp]
            } else {
              beta_temp[[j]][labels_temp[[j]][i, t], t] <- newclustervalue
              sort_beta_temp <- unique(labels_temp[[j]][, t])
              beta_new <- beta_temp[[j]][sort_beta_temp, t]
            }
          } else {
            sort_beta_temp <- unique(labels_temp[[j]][, t])
            beta_new <- beta_temp[[j]][sort_beta_temp, t]
          }
          
          beta_temp[[j]][, t] <- NA
          beta_temp[[j]][1:length(sort_beta_temp), t] <- beta_new
          labels_temp[[j]][, t] <- labels_new
        }
        
        # --- Update Spline Coefficients (Beta) ---
        n_cluster <- max(labels_temp[[j]][, t])
        Y_t <- t(vapply(Y, function(y) y[t, ], FUN.VALUE = double(length(ages))))
        beta_t <- vapply(beta_temp, function(b) b[, t], FUN.VALUE = double(n))
        labels_t <- vapply(labels_temp, function(lab) lab[, t], FUN.VALUE = double(n))
        
        for (k in 1:n_cluster){
          beta_temp[[j]][k, t] <- up_beta(j = j, 
                                          k = k, 
                                          t = t,
                                          Y_t = Y_t,
                                          beta_t = beta_t,
                                          phi_jt = phi_temp[[j]][t], 
                                          delta_j = delta_temp[[j]],
                                          sigma_vec = sigma_temp,
                                          labels_t = labels_t,
                                          spline_basis = S)
          beta_t[k] <- beta_temp[[j]][k, t]
        }
        
      } # End loop over t
      
      # Store parameter states for iteration d
      labels_res[[j]][,, d] <- as.integer(labels_temp[[j]])
      gamma_res[[j]][,, d] <- gamma_temp[[j]]
      beta_res[[j]][,, d] <- beta_temp[[j]]
      
      # --- Update Delta ---
      to_up_deltaj <- (t(beta_temp[[j]]) - GP_means[[j]])[!is.na(t(beta_temp[[j]]))]
      delta_res[[j]][d] <- delta_temp[[j]] <- 
        sqrt(1 / rgamma(1, 
                        shape = a_delta + 0.5 * length(to_up_deltaj),
                        rate = b_delta + 0.5 * sum(to_up_deltaj^2)))
      
      # --- Update Phi ---
      prec.lik_phi <- diag(apply(beta_temp[[j]], 2, 
                                 function(x) sum(!is.na(x))) / delta_temp[[j]]^2, 
                           nrow = T_final)
      var.post_phi <- mysolve(SIGMAinv / omega_temp[[j]]^2 + prec.lik_phi)
      phi_res[[j]][, d] <- phi_temp[[j]] <- 
        drop(rmvn(1,
                  mu = var.post_phi %*% ((SIGMAinv / omega_temp[[j]]^2) %*% GP_means[[j]] + 
                                           colSums(beta_temp[[j]] / delta_temp[[j]]^2, na.rm = TRUE)),
                  sigma = var.post_phi))
      
      # --- Update Omega ---
      omega_res[[j]][d] <- omega_temp[[j]] <- 
        sqrt(1 / rgamma(1, 
                        shape = a_omega + 0.5 * T_final,
                        rate = b_omega + 0.5 * sum((t(cholSIGinv) %*% (phi_temp[[j]] - GP_means[[j]]))^2)))
      
      # --- Update Alpha ---
      alpha_res[[j]][d] <- alpha_temp[[j]] <- 
        up_alpha_j(sum(gamma_temp[[j]]), n * (T_final -1), a_alpha, b_alpha)
      
      # --- Update CRP Concentration Parameter (M) ---
      M_res[[j]][d] <- M_temp[[j]] <- up_M(gamma_j = gamma_temp[[j]],
                                           labels_j = labels_temp[[j]],
                                           M_j = M_temp[[j]],
                                           a_M = a_M, b_M = b_M,
                                           T_final = T_final, n = n)
      
    } # End loop over j
    
    # --- Update Observation Variances (Sigma) ---
    means <- vapply(1:T_final,
                    function(t) vapply(1:p, 
                                       function(x) beta_temp[[x]][labels_temp[[x]][, t], t],
                                       FUN.VALUE = double(n)) %*% t(S),
                    FUN.VALUE = matrix(1, n, length(ages)))
    
    for (i in 1:n){
      sigma_res[i, d] <- sigma_temp[i] <- sqrt(1 / rgamma(1, 
                                                          shape = a_sigma + 0.5 * T_final * length(ages),
                                                          rate = b_sigma + 0.5 * sum((Y[[i]] - t(means[i,, ])) ^ 2)))
    }
    
    # Verbose logging
    if ((d %% print_step) == 1) {
      now <- Sys.time()
      cat("iter. ", d - 1, " completed",
          "\n Time elapsed since start: ", difftime(now, inizio, units = "mins"),
          "\n Time elapsed since last log: ", difftime(now, last_time, units = "mins"), " mins", 
          "\n Current time: ", format(Sys.time(), "%H:%M:%S"), 
          "\n\n", sep = "")
      last_time <- now
    }
    
  } # End loop over d
  
  fine <- Sys.time()
  exec_time <- difftime(fine, inizio)
  
  out <- list("res" = list("labels" = labels_res,
                           "beta" = beta_res,
                           "delta" = delta_res,
                           "phi" = phi_res,
                           "omega" = omega_res,
                           "alpha" = alpha_res,
                           "gamma" = gamma_res,
                           "M" = M_res,
                           "sigma" = sigma_res,
                           "exec_time" = exec_time),
              "parameters" = list("a_delta" = a_delta,
                                  "b_delta" = b_delta,
                                  "a_omega" = a_omega,
                                  "b_omega" = b_omega,
                                  "a_alpha" = a_alpha,
                                  "b_alpha" = b_alpha,
                                  "a_M" = a_M,
                                  "b_M" = b_M,
                                  "a_sigma" = a_sigma,
                                  "b_sigma" = b_sigma),
              "data" = Y,
              "other" = list("n" = n,
                             "T_final" = T_final,
                             "p" = p,
                             "n_iter" = n_iter,
                             "ages" = ages,
                             "S" = S,
                             "path_save" = path_save,
                             "name_save" = name_save))
  
  saveRDS(out, paste0(path_save, name_save))
  return(out)
}
