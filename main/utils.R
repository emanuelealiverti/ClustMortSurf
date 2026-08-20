############################################################################
## Utility functions for Gibbs sampler and results extraction of          ##
## "Bayesian local clustering of age-period mortality surfaces            ##
## across multiple countries" (Romanò, Aliverti and Durante)              ##
############################################################################

## ----------------------------------------------------------------------
## Load all required libraries
## ----------------------------------------------------------------------
require(ggplot2)
require(mvnfast)
require(mcclust)
require(mcclust.ext)
require(dplyr)
require(purrr)
require(tidyr)
require(reshape2)
require(multvardiv)
require(aricode)


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
## Computes the posterior similarity matrix (PSM) for each spline basis
## and time point, measuring the pairwise co-clustering probabilities across
## MCMC samples.
##
## Arguments
##    import    logical; if TRUE, loads pre-computed PSM from "psm.RDS"
##    save      logical; if TRUE, saves the computed PSM to "psm.RDS"
##    draws     list of MCMC sample arrays for cluster allocations
##    burnin    number of initial MCMC iterations to discard as burn-in
##
## Value: a list of 3D arrays containing co-clustering probabilities for each spline basis.

compute_prob_coclust <- function(import = FALSE, save = FALSE, draws = NULL, burnin = NULL){


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
## Extracts a point estimate of the partition structure by minimizing the
## Variation of Information (VI) loss function over MCMC draws given a PSM.
##
## Arguments
##    import    logical; if TRUE, loads pre-computed point estimates from "ppe.RDS"
##    save      logical; if TRUE, saves the computed point estimates to "ppe.RDS"
##    draws     list of MCMC sample arrays for cluster allocations
##    psm       list of posterior similarity matrices obtained from compute_prob_coclust
##    burnin    number of initial MCMC iterations to discard as burn-in
##
## Value: a list of matrices with cluster label point estimates per unit and time point.

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



## ----------------------------------------------------------------------
## Minimize Variation of Information with multiple restarts
## ----------------------------------------------------------------------
## Helper function that runs VI minimization using 'mcclust.ext' initialized
## at multiple random partition states to avoid local minima.
##
## Arguments
##    draws    matrix of posterior cluster allocation draws
##    psm      posterior similarity matrix; if NULL, computed on the fly
##
## Value: a vector of optimal cluster labels minimizing VI.
my.minVI <- function(draws, psm = NULL){

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
## Evaluates co-clustering accuracy against a ground truth partition. At each
## MCMC iteration, the sampled partition is compared to the true one
## for every spline basis and time point, obtaining a p x T matrix with the
## proportion of correctly co-clustered pairs of units. The final accuracy is
## the average of these matrices over the MCMC sample.
## 
##
## Arguments
##    draws            list of MCMC sample arrays for cluster allocations (units x years x MCMC iterations) (discarding burnin)
##    clusters_true    list of matrices containing true cluster memberships
##    burnin           number of initial MCMC iterations to discard as burn-in
##
## Value: a list containing the matrix of posterior mean co-clustering accuracy ('psm_acc').

summarise_clusters <- function(draws_full, clusters_true, burnin){


	draws  <- lapply(draws_full, \(.l) .l[,,-burnin])
	n <- dim(draws[[1]])[1]
	TT <- dim(draws[[1]])[2]
	p <- length(draws)
	iters <- seq_len(dim(draws[[1]])[3])

	# True co-clustering adjacency
	# List of p matrices of size n x n x TT
	adj_true <- vector("list", p)
	for (j in seq_len(p)){
		adj_true[[j]] <- vapply(seq_len(TT), function(t){
			a <- 1*outer(clusters_true[[j]][ , t], clusters_true[[j]][ , t], "==")
			diag(a) <- 0
			a
		}, matrix(0, n, n))
	}

	psm_acc <- matrix(0, p, TT) # computed recursively
	dimnames(psm_acc) <- list("Spline" = paste0("Spline", seq_len(p)),
				  "Year" = 1:TT)

	for (it in iters){

		# Proportion of correctly co-clustered pairs at this iteration, for every (spline, year)
		acc_it <- matrix(NA, p, TT)
		for (j in seq_len(p)){
			labs_it <- draws[[j]][ , , it]
			for (t in seq_len(TT)){
				adj_it <- 1*outer(labs_it[ , t], labs_it[ , t], "==")
				coclust_it  <- (adj_it == adj_true[[j]][,,t])
				acc_it[j, t] <- mean(coclust_it[lower.tri(coclust_it)])
			}
		}

		psm_acc <- psm_acc + acc_it
	}

	psm_acc <- psm_acc / length(iters)

	return(psm_acc)

}



## ----------------------------------------------------------------------
## Obtain estimates for individual spline coefficients from MCMC samples
## ----------------------------------------------------------------------
## Maps cluster-level spline coefficients back to individual observation units
## by matching cluster draws with unit assignment draws across MCMC iterations.
##
## Arguments
##    import        logical; if TRUE, loads pre-computed coefficients from "beta_units.RDS"
##    save          logical; if TRUE, saves the computed unit coefficients to file
##    draws_beta    list of MCMC sample arrays for cluster-level spline coefficients
##    draws_labs    list of MCMC sample arrays for unit cluster assignments
##    burnin        number of initial MCMC iterations to discard as burn-in
##
## Value: a list of three 3D arrays (units x splines x years) with the posterior
## means ('Post_Mean') and the lower ('Post_q025') and upper ('Post_q975')
## endpoints of the 95% equal-tailed credible interval of the unit-level spline
## coefficients.

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

		empty_array <- array(NA, dim = c(n, p, TT),
				     dimnames = list(paste0("Unit", seq_len(n)),
						     paste0("Spline", seq_len(p)),
						     seq_len(TT)))
		beta_units <- list(Post_Mean = empty_array, Post_q025 = empty_array, Post_q975 = empty_array)

		for (j in 1:p){

			cat(j, "\n")

			beta <- draws_beta[[j]][ , , -burnin]
			labs <- draws_labs[[j]][ , , -burnin]


			for (t in seq_len(TT)){
				l <- labs[ , t, ]
				b <- beta[ , t, ]
				# Work with vector version of b and l for faster execution
				l_vec = c(l) + rep((nrow(l))*(1:ncol(l)-1), each = nrow(l))
				b_vec = c(b)
				# Unit-level draws of beta_ijt: units on rows, retained iterations on columns
				draws_units <- matrix(b_vec[l_vec], nrow = nrow(l), ncol = ncol(l))

				beta_units$Post_Mean[, j, t] <- rowMeans(draws_units)
				ci <- apply(draws_units, 1, quantile, probs = c(0.025, 0.975), names = FALSE)
				beta_units$Post_q025[, j, t] <- ci[1, ]
				beta_units$Post_q975[, j, t] <- ci[2, ]
			}

		}
		if (save){saveRDS(beta_units, "beta_units.RDS")}
	}

	return(beta_units)
}



## ----------------------------------------------------------------------
## Summarise estimates of individual spline coefficients
## ----------------------------------------------------------------------
## Merges estimated spline coefficients with true benchmark coefficients into
## a single tidy data frame for post-estimation comparison and plotting.
##
## Arguments
##    estimated    list with the 3D arrays 'Post_Mean', 'Post_q025' and 'Post_q975'
##                 of estimated unit spline coefficients (from compute_beta_units)
##    true         3D array of true unit spline coefficients
##
## Value: a data frame joining estimated and true coefficients by Unit, Spline and
## Year, with the posterior mean ('Post_Mean'), the 95% credible interval
## endpoints ('Post_q025', 'Post_q975') and the true value ('Beta_true').

summarise_beta <- function(estimated, true){


	stopifnot(all(c("Post_Mean", "Post_q025", "Post_q975") %in% names(estimated)))

	names(dimnames(true)) <- c("Unit", "Spline", "Year")

	est_df <- map_dfr(estimated, ~ as.data.frame.table(.x), .id = "index") %>%
		pivot_wider(names_from = index, values_from = Freq) %>%
		rename(Unit = Var1, Spline = Var2, Year = Var3) %>%
		mutate(Year = as.numeric(Year))

	beta_df <- full_join(est_df,
			     melt(true, value.name = "Beta_true"),
			     by = c("Unit", "Spline", "Year"))

	return(beta_df)
}
