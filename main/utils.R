############################################################################
## Utility functions for Gibbs sampler and results extraction of          ##
## "Bayesian local clustering of age-period mortality surfaces            ##
## across multiple countries" (Romanò, Aliverti and Durante)              ##
############################################################################

require(mvnfast)



## ----------------------------------------------------------------------
## Invert PD matrices (e.g. varcov)
## ----------------------------------------------------------------------
##
## Arguments
##    x   PD matrix
##
## Value: inverted matrix.

mysolve <- function(x){
  
  inv_x = chol2inv(chol(x))
  
  return(inv_x)
}



## ----------------------------------------------------------------------
## Partition representations conversion: labels -> sets
## ----------------------------------------------------------------------
##
## Arguments
##    labels    vector of membership labels
##
## Value: list of vectors corresponding to sets in partition.

lab2rho <- function(labels){
  
  # Number of tables 
  nTbl <- max(labels)
  
  # Partition sets representation
  rho <- lapply(1:nTbl, function(t) which(labels == t)) 
  
  return(rho)
}



## ----------------------------------------------------------------------
## Partition representations conversion: sets -> labels
## ----------------------------------------------------------------------
##
## Arguments
##    rho   list of vectors corresponding to sets in partition 
##
## Value: vector of membership labels.

rho2lab <- function(rho){
  
  # Number of entities in the partition
  n <- max(unlist(rho))
  
  # Partition labels representation
  labels <- vapply(1:n, 
                   function(x) which(vapply(rho, is.element, el = x, FUN.VALUE = FALSE)), 
                   FUN.VALUE = integer(1))
  
  return(labels)
}



## ----------------------------------------------------------------------
## CRP pmf
## ----------------------------------------------------------------------
## Arguments
##    size    vector with cluster sizes
##    M       CRP concentration parameter
##
## Value: probability of observing a partition with cluster sizes encoded in the
## vector size.

dCRP <- function(size, M){
  
  require(mggd)
  
  # Remove empty clusters, if any
  size <- size[size>0]
  
  # Compute CRP pmf
  pmf = M^NROW(size) * prod(factorial(size - 1)) / pochhammer(M, sum(size))
  
  return(pmf)
}



## ----------------------------------------------------------------------
## CRP sampler
## ----------------------------------------------------------------------
##
## Arguments
##    n   sample size
##    M   CRP concentration parameter
##
## Value: a list with the sets representation of the sampled partition.
rCRP <- function(n, M){
  
  # Initialize labels vector
  labels <- rep(0, n)
  
  # Allocation first unit
  labels[1] <- 1
  
  for (dnr in 2:n) {
    
    # Compute occupation probabilities for current diner
    vOcc <- table(labels[1:(dnr-1)])
    vProb <- c(vOcc, M) / (dnr - 1 + M)
    
    # Add table label to diner
    nTbl <- as.numeric(names(vOcc)[length(vOcc)]) 
    labels[dnr] <- sample.int(nTbl+1, size=1, prob=vProb)
  }
  
  nTbl <- max(c(nTbl, labels[n]))
  rho <- lapply(1:nTbl, function(t) which(labels == t))
  
  return(rho)
}



## ----------------------------------------------------------------------
## Squared exponential kernel
## ----------------------------------------------------------------------
## Evaluation of the squared exponential kernel used in the GP varcov matrix.
## The kernel is computed at distance d between two time points t1 and t2.
## The term to control the scale (which is included for example in BDA3) is 
## not considered here because it is considered outside the kernel and a prior 
## is placed on it.
## 
## Arguments
##
##    l   smoothness parameter (the larger l and the greater the smoothness)
##    d   |t1 - t2|
## Value: a scalar with the evalution of the squared exponential at d.

sq_exp_ker <- function(d, l){
  
  sq_eval = exp(- d^2 / (2*l^2))
  
  return(sq_eval)
}


## ----------------------------------------------------------------------
## Compute posterior similarity matrix from MCMC draws
## ----------------------------------------------------------------------
## DESCRIPTION TO BE DONE.
##
## Arguments
##
##    import
##    save
##    draws
##    burnin
##
## Value: .
compute_prob_coclust <- function(import = FALSE, save = FALSE, draws = NULL, burnin = NULL){
  
  require(mcclust)
  
  if (!import & is.null(draws)){stop("Please provide MCMC draws or select import = TRUE.")}
  if (import & !is.null(draws)){warning("MCMC draws are provided even if import = TRUE.") }
  
  if (import){
    
    prob_coclust <- readRDS("psm.RDS")
    
  } else {
    
    n <- nrow(draws[[1]])
    TT <- ncol(draws[[1]])
    p <- length(draws)
    
    prob_coclust <- vector("list", p)
    names(prob_coclust) <- paste0("Spline", seq_len(p))
    
    for (j in 1:p){
      cat(j, "\t")
      
      ppc <- vapply(1:TT, 
                    function(t) comp.psm(t(draws[[j]][, t, -burnin])),
                    matrix(1, n, n))
      
      dimnames(ppc) <- list("Units" = paste0("Unit", seq_len(n)),
                            "Units" = paste0("Unit", seq_len(n)),
                            "Year" = 1:TT)
      
      prob_coclust[[j]] <- ppc
    }
    
    if (save){saveRDS(prob_coclust, "psm.RDS")}
  }
  
  return(prob_coclust)
}



## ----------------------------------------------------------------------
## Compute partition point estimate from PSM
## ----------------------------------------------------------------------
## DESCRIPTION TO BE DONE.
##
## Arguments
##
##    import
##    save
##    draws
##    psm
##    burnin
##
## Value: .

compute_partition_point_est <- function(import = FALSE, save = FALSE, 
                                        draws = NULL, psm = NULL, burnin = NULL){
  
  if (!import & is.null(draws)){stop("Please provide MCMC draws or select import = TRUE.")}
  if (import & (!is.null(draws) | !is.null(psm))){
    warning("MCMC draws and/or PSM are provided even if import = TRUE.")
  }
  
  if (import){
    
    partition_point_est <- readRDS("ppe.RDS")
    
  } else {
    
    n <- dim(draws[[1]])[1]
    TT <- dim(draws[[1]])[2]
    p <- length(draws)
    
    partition_point_est <- vector("list", p)
    names(partition_point_est) <- paste0("Spline", seq_len(p))
    
    for (j in 1:p){
      cat(j, "\t")
      
      ppe <- vapply(1:TT, 
                    function(t) my.minVI(t(draws[[j]][, t, -burnin]), psm[[j]][ , , t]),
                    integer(n))
      
      dimnames(ppe) <- list("Units" = paste0("Unit", seq_len(n)),
                            "Year" = 1:TT)
      
      partition_point_est[[j]] <- ppe
    }
    
    if (save){saveRDS(partition_point_est, "ppe.RDS")}
  }
  
  return(partition_point_est)
}


my.minVI <- function(draws, psm = NULL){
  
  require(mcclust)
  require(mcclust.ext)
  
  if (is.null(psm)){psm <- comp.psm(draws)}
  start_set <- draws[sample(1:nrow(draws), 10, replace = FALSE), ]
  
  list_opt <- apply(start_set,
                    1,
                    function(start) minVI(psm = psm, 
                                          method = "avg", 
                                          max.k = nrow(psm),
                                          start.cl = start))
  pick <- which.max(sapply(list_opt, function(x) x$value))
  out <- list_opt[[pick]]$cl
  
  return(out)
}


## ----------------------------------------------------------------------
## Summarise estimates of cluster labels
## ----------------------------------------------------------------------
## DESCRIPTION TO BE DONE.
##
## Arguments
##
##    ppe
##    psm
##    clusters_true
##
## Value: .

summarise_clusters <- function(ppe, psm, clusters_true){
  
  require(aricode)
  
  n <- nrow(ppe[[1]])
  TT <- ncol(ppe[[1]])
  p <- length(ppe)
  
  nvi <- matrix(NA_real_, p, TT)
  psm_acc <- matrix(NA_real_, p, TT)
  dimnames(nvi) <- dimnames(psm_acc) <- list("Spline" = paste0("Spline", seq_len(p)),
                                             "Year" = 1:TT)
  
  for (j in seq_len(p)){
    for (t in seq_len(TT)){
      nvi[j, t] <- NVI(clusters_true[[j]][ , t], ppe[[j]][ , t])
      adj_true <- 1*outer(clusters_true[[j]][ , t], clusters_true[[j]][ , t], "==")
      psm_acc[j, t] <- mean(adj_true*psm[[j]][ , , t] + (1-adj_true)*(1-psm[[j]][ , , t]))
    }
  }
  
  return(list(psm_acc = psm_acc, ppe_acc = 1 - nvi))
  
}



## ----------------------------------------------------------------------
## Obtain estimates for individual spline coefficients from MCMC samples
## ----------------------------------------------------------------------
## DESCRIPTION TO BE DONE.
##
## Arguments
##
##    import
##    save
##    draws_beta
##    draws_labs
##    burnin
##
## Value: .
compute_beta_units <- function(import = FALSE, save = FALSE, 
                               draws_beta = NULL, draws_labs = NULL, burnin = NULL){
  
  
  if (!import & (is.null(draws_beta) | is.null(draws_labs))){
    stop("Please provide MCMC draws or select import = TRUE.")
  }
  if (import & (!is.null(draws_labs) | !is.null(draws_beta))){
    warning("MCMC draws are provided even if import = TRUE.")
  }
  
  if (import){
    
    beta_units <- readRDS("beta_units.RDS")
    
  } else {
    
    n <- dim(draws_labs[[1]])[1]
    TT <- dim(draws_labs[[1]])[2]
    p <- length(draws_labs)
    
    beta_units  <- array(NA_real_, dim = c(n, p, TT),
                         dimnames = list(paste0("Unit", seq_len(n)),
                                         paste0("Spline", seq_len(p)),
                                         seq_len(TT)))
    
    for (j in 1:p){
      
      cat(j, "\t")
      
      beta <- draws_beta[[j]][ , , -burnin]
      labs <- draws_labs[[j]][ , , -burnin]
      
      
      for (t in seq_len(TT)){
          l <- labs[ , t, ]
          b <- beta[ , t, ]
          # Work with vector version of b and l for faster execution
          l_vec = c(l) + rep((nrow(l))*(1:ncol(l)-1), each = nrow(l))
          b_vec = c(b)
          beta_units[, j, t] <- rowMeans(matrix(b_vec[l_vec], nrow = nrow(l), ncol = ncol(l)))
      }
      
    }
    if (save){saveRDS(beta_units, paste0(gender, "/beta_units.rds"))}
  }
  
  return(beta_units)
}



## ----------------------------------------------------------------------
## Summarise estimates of individual spline coefficients
## ----------------------------------------------------------------------
## DESCRIPTION TO BE DONE.
##
## Arguments
##
##    ppe
##    psm
##    clusters_true
##
## Value: .
summarise_beta <- function(estimated, true){
  
  names(dimnames(estimated)) = names(dimnames(true)) = c("Unit", "Spline", "Year")
  beta_df <- full_join(melt(estimated, value.name = "Beta_est"),
                       melt(true, value.name = "Beta_true"),
                       by = c("Unit", "Spline", "Year"))
  
  
}