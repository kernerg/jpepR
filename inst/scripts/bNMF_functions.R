### FUNCTIONS

# ******************************************************
# Udler et al. bNMF MODEL (PLEIOTROPIC VERSION 0)
# ******************************************************

BayesNMF.L2EU <- function(
    V0, n.iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0 #20, 10
) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits)
  # n.iter: Number of iterations for parameter optimization
  # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
  # phi: Scaling parameter
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  W <- matrix(runif(N * K) * Vmax, ncol=K)
  H <- matrix(runif(M * K) * Vmax, ncol=M)
  
  I <- array(1, dim=c(N, M))
  V.ap <- W %*% H + eps
  
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    # H <- H * (t(W) %*% V) / 
    #   (t(W) %*% V.ap + phi * H * matrix(rep(1 / rep(1,K), M), ncol=M) + eps)
    # V.ap <- W %*% H + eps
    # W <- W * (V %*% t(H)) / 
    #   (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / rep(1,K), N), ncol=N)) + eps)
    
    H <- H * (t(W) %*% V) /
      (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol=M) + eps)
    V.ap <- W %*% H + eps
    W <- W * (V %*% t(H)) /
      (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol=N)) + eps)
    
    V.ap <- W %*% H + eps
    
    # # Introduce regularization for elements outside the desired range
    # V.ap[V.ap < 1e-4] <- 0
    # V.ap[V.ap > 1] <- 1
    
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    #del <- sum((V - V.ap)^2) / 2
    like <- sum((V - V.ap)^2) / 2
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / 
                                         lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del, 
          sum(colSums(W) != 0), sum(lambda >= lambda.cut), '\n')
    }
    iter <- iter + 1
  }
  return(list(
    W,  # Variant weight matrix (N x K)
    H,  # Trait weight matrix (K x M)
    n.like,  # List of reconstruction errors (sum of squared errors / 2) per iteration
    n.evid,  # List of negative log-likelihoods per iteration
    n.lambda,  # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
    n.error  # List of reconstruction errors (sum of squared errors) per iteration
  ))
}


run_bNMF <- function(z_mat, n_reps=10, random_seed=1, K=20, K0=10, tolerance=1e-7) {
  
  # Given an input matrix as created by prep_z_matrix(), run the bNMF procedure
  # a series of times to generate results and evaluate cluster stability
  
  print(paste0("Running bNMF clustering procedure (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!",tolerance))
  
  set.seed(random_seed)
  
  bnmf_reps <- lapply(1:n_reps, function(r) {
    print(paste("ITERATION",r))
    res <- BayesNMF.L2EU(V0 = z_mat, K=K, K0=K0, tol=tolerance)
    names(res) <- c("W", "H", "n.like", "n.evid", "n.lambda", "n.error")
    res
  })
  bnmf_reps
}


make_run_summary <- function(reps) {
  
  # Given a list of bNMF iteration outputs, summarize the K choices and associated likelihoods across runs
  
  run_summary <- map_dfr(1:length(reps), function(i) {
    res <- reps[[i]]
    final_lambdas <- res$n.lambda[[length(res$n.lambda)]]
    tibble(
      run=i,
      K=sum(final_lambdas > min(final_lambdas)),  # Assume that lambdas equal to the minimum lambda are ~ 0
      evid=res$n.evid[[length(res$n.evid)]]  # Evidence = -log_likelihood
    )
  }) %>%
    arrange(evid)
  
  unique.K <- table(run_summary$K)
  n.K <- length(unique.K)  # Number of distinct K
  MAP.K.run <- sapply(names(unique.K), function(k) {  # bNMF run index with the maximum posterior for given K
    tmp <- run_summary[run_summary$K == k, ]
    tmp$run[which.min(tmp$evid)]
  })
  
  list(run_tbl=run_summary, unique.K=unique.K, MAP.K.run=MAP.K.run)
}

# ******************************************************
# bNMF PLEIOTROPIC
# ******************************************************

BayesNMF.L2EU.trait <- function(
    V0, W_init = NULL, sp = FALSE, n.iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0 #20, 10
) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits)
  # n.iter: Number of iterations for parameter optimization
  # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
  # phi: Scaling parameter
  # sp: indicates whether or not to use sparsity penalty
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  if (is.null(W_init)) {
    W <- matrix(runif(N * K) * Vmax, ncol=K)
  } else {
    W <- W_init
  }
  H <- matrix(runif(M * K) * Vmax, ncol=M)
  
  I <- array(1, dim=c(N, M))
  V.ap <- W %*% H + eps
  
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% V) /
      (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol=M) + eps)
    
    # enforces sparsity by associating each trait to only one cluster
  #  if(sp) {
  #    for (i in 1:ncol(H)) {
  #      max_val_index <- which.max(H[, i])
  #      H[-max_val_index, i] <- 0  # Set all other entries in the row to 0
  #    }
  #  }
    
    V.ap <- W %*% H + eps
    W <- W * (V %*% t(H)) /
      (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol=N)) + eps)
    
    V.ap <- W %*% H + eps
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    like <- sum((V - V.ap)^2) / 2
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / 
                                         lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del, 
          sum(colSums(W) != 0), sum(lambda >= lambda.cut), '\n')
    }
    iter <- iter + 1
  }

  return(list(
    W = W,
    H = H,
    n.like = n.like,
    n.evid = n.evid,
    n.lambda = n.lambda,
    n.error = n.error
  ))
}

run_BayesNMF.trait <- function(V0, W_init = NULL, sp = FALSE, n.iter = 2000,
                         a0 = 10,
                         tol = 1e-7, K = 10, K0 = 10,
                         phi = 1.0, n_reps,
                         random_seed = 1) {
 
  print(paste0("Running bNMF clustering procedure (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!",tol))
    
  set.seed(random_seed)
  successful_runs <- list()
  while (length(successful_runs) < n_reps) {
    try({
      print(paste("Attempt:", length(successful_runs) + 1))
      res <- BayesNMF.L2EU.trait(
        V0 = V0, W_init = W_init, sp = sp, n.iter = n.iter, a0 = a0,
        tol = tol, K = K, K0 = K0, phi = phi)
      successful_runs[[length(successful_runs) + 1]] <- res
    }, silent = TRUE)
  }
  return(successful_runs)
}


# ******************************************************
# bNMF EPIGENOMIC
# ******************************************************


BayesNMF.L2EU.tissue <- function(
    V0, W_init = NULL, sp = FALSE, n.iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0 #20, 10
) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits)
  # n.iter: Number of iterations for parameter optimization
  # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
  # phi: Scaling parameter
  # sp: indicates whether or not to use sparsity penalty
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  if (is.null(W_init)) {
    W <- matrix(runif(N * K) * Vmax, ncol=K)
  } else {
    W <- W_init
  }
  H <- matrix(runif(M * K) * Vmax, ncol=M)
  
  I <- array(1, dim=c(N, M))
  V.ap <- W %*% H + eps
  
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% V) /
      (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol=M) + eps)
    
    # enforces sparsity by associating each trait to only one cluster
   if(sp) {
     for (i in 1:ncol(H)) {
       max_val_index <- which.max(H[, i])
       H[-max_val_index, i] <- 0  # Set all other entries in the row to 0
     }
   }
    
    V.ap <- W %*% H + eps
    W <- W * (V %*% t(H)) /
      (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol=N)) + eps)
    
    V.ap <- W %*% H + eps
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    like <- sum((V - V.ap)^2) / 2
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / 
                                         lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del, 
          sum(colSums(W) != 0), sum(lambda >= lambda.cut), '\n')
    }
    iter <- iter + 1
  }

  return(list(
    W = W,
    H = H,
    n.like = n.like,
    n.evid = n.evid,
    n.lambda = n.lambda,
    n.error = n.error
  ))
}

run_BayesNMF.tissue <- function(V0, W_init = NULL, sp = FALSE, n.iter = 2000,
                         a0 = 10,
                         tol = 1e-7, K = 10, K0 = 10,
                         phi = 1.0, n_reps,
                         random_seed = 1) {
 
  print(paste0("Running bNMF clustering procedure (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!",tol))
    
  set.seed(random_seed)
  successful_runs <- list()
  while (length(successful_runs) < n_reps) {
    try({
      print(paste("Attempt:", length(successful_runs) + 1))
      res <- BayesNMF.L2EU.tissue(
        V0 = V0, W_init = W_init, sp = sp, n.iter = n.iter, a0 = a0,
        tol = tol, K = K, K0 = K0, phi = phi)
      successful_runs[[length(successful_runs) + 1]] <- res
    }, silent = TRUE)
  }
  return(successful_runs)
}

# ******************************************************
# bNMF JOINT (aka JPEP), VERSION 1
# ******************************************************

BayesNMF.L2EU.joint <- function(
    V0, V0_tis, W_init = NULL, sp = FALSE, alpha_proj, alpha_cross,
    n.iter = 2000, tol = 1e-7, K = 15, alpha_tis = 0.5
) {

  # Bayesian NMF with half-normal priors for W and H
  # V0: input pip matrix (variants x traits)
  # Vtis: input tissue matrix (variants x tissues)
  # n.iter: Number of iterations for parameter optimization
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)

  # Initialize variables
  eps <- 1.e-50; del <- 1.0
  V <- V0 - min(V0)
  V_tis <- V0_tis - min(V0_tis)
  N <- nrow(V); M <- ncol(V); M_tis <- ncol(V_tis)
  Vmax <- max(V); V_tismax <- max(V_tis)

  # Initialize hidden factors
  W <- if (is.null(W_init)) matrix(runif(N * K) * Vmax, ncol = K) else W_init
  H <- matrix(runif(M * K) * Vmax, ncol = M)
  H_tis <- matrix(runif(M_tis * K) * V_tismax, ncol = M_tis)

  # Initialize reconstructed input matrices
  V.ap <- W %*% H + eps
  V.ap_tis <- W %*% H_tis + eps

  # Initialize penalization terms
  param_non_penalized_clusters <- 1e10
  worst_cluster <- 0.99
  cluster_to_not_touch <- 5e-2
  lambda_proj <- rep(param_non_penalized_clusters, nrow(H))
  lambda_cross <- rep(param_non_penalized_clusters, nrow(H))

  # Initialize projection matrices
  V_proj <- V %*% t(H)
  V_tis_proj <- V_tis %*% t(H_tis)

  # intialize output lists
  n.lambda_proj <- list(); n.lambda_proj[[1]] <- lambda_proj
  n.lambda_cross <- list(); n.lambda_cross[[1]] <- lambda_cross

  iter <- 2; p_cross <- NA
  while ((del >= tol) &&  iter < n.iter) {

    # Update matrices
    H <- update_H(W, V, V.ap, H, lambda_proj, lambda_cross, alpha_proj, alpha_cross)
    H_tis <- update_H_tis(W, V_tis, V.ap_tis, H_tis, lambda_proj, lambda_cross, alpha_proj, alpha_cross)
    V.ap <- W %*% H + eps
    W <- update_W_joint(W, V, V.ap, V_tis, V.ap_tis, H, H_tis,
                        lambda_proj, lambda_cross, alpha_tis, alpha_proj, alpha_cross)

    # Update reconstructed matrices
    V.ap <- W %*% H + eps
    V.ap_tis <- W %*% H_tis + eps

    # Current surviving clusters
    myk <- which(rowSums(H) > 1e-10 & rowSums(H_tis) > 1e-10)
    if (length(myk) < 2) { stop("Less than 2 clusters") }

    # Update projection matrices
    if (alpha_proj > 0 || alpha_cross > 0) {
      V_proj <- V %*% t(H)
      V_tis_proj <- V_tis %*% t(H_tis)
    }

    # Compute correlation between projections at each cluster (p values are used for penalty)
    # And penalize clusters with low projection correlation
    if (alpha_proj > 0) {
      cors <- calculate_correlations(V_proj, V_tis_proj, H, H_tis)
      correlations <- cors$correlations
      ps <- cors$p
    }

    # Compute cross-correlations between projections (p values are used for penalty)
    # And penalize clusters that have strong projection correlation with another cluster
    #if (alpha_cross > 0) {
      cross_mat <- pmin(calculate_cross_correlations(V_proj,V_tis_proj, H, H_tis) /
        max(correlations), 1)
      if (any(is.na(cross_mat))) { stop("Less than 2 clusters") }
      cross_matp <- calculate_cross_correlationsp(V_proj,V_tis_proj, H, H_tis)
      p_cross <- apply(cross_matp, 1, min)

      # Define groups of clusters to eventually impose penalty on all cluster except one of each group
      graph <- graph_from_adjacency_matrix(abs(cross_mat), mode = "undirected", weighted = TRUE)
      # Apply the Louvain community detection algorithm
      communities <- cluster_louvain(graph)
      community_groups <- communities$membership
      g <- as.numeric(names(table(community_groups)[table(community_groups) > 1]))
      # Retain the cluster with highest tissue+trait variance at each defined group of clusters
      v <- apply(H, 1, var)
      v_tis <- apply(H_tis, 1, var)
      # Apply penalty to the rest of clusters
      punish <- c(); punish_myk <- c()
      for (gr in g) {
        cl <- myk[which(community_groups == gr)]
       p <- cl[which.max(v[cl] + v_tis[cl])]
        punish <- c(punish, setdiff(cl, p))
       punish_myk <- c(punish_myk, which(myk %in% setdiff(cl, p)))
      }
    #}

    # Impose penalty on selected clusters
    # alreadygood <- which(ps < cluster_to_not_touch)
    # lambda_proj[alreadygood] <- param_non_penalized_clusters
    if (alpha_proj > 0) {
      lambda_proj <- rep(-log10(worst_cluster), nrow(H))
      lambda_proj[myk] <- (-log10(pmin(ps, worst_cluster)))
    }
    
    if (alpha_cross > 0) {
      lambda_cross <- rep(param_non_penalized_clusters, nrow(H))
      lambda_cross[punish] <- 1 / (-log10(p_cross[punish_myk]))
    }

    if (length(lambda_proj) == length(n.lambda_proj[[iter - 1]])) {
      del1 <- max(abs(lambda_proj - n.lambda_proj[[iter - 1]]) / n.lambda_proj[[iter - 1]])
      del2 <- max(abs(lambda_cross - n.lambda_cross[[iter - 1]]) / n.lambda_cross[[iter - 1]])
      del <- ifelse(alpha_cross > 0, max(del1, del2), del1)
    }

    # save penalties to define convergence
    n.lambda_proj[[iter]] <- lambda_proj
    n.lambda_cross[[iter]] <- lambda_cross

    if (iter %% 100 == 0) {
      cat(iter, del, sum(colSums(W) != 0), "\n")
    }
    iter <- iter + 1
  }

  if (iter >= n.iter) {
    stop("Maximum number of iterations reached without convergence.")
  }

  return(list(
    W = W, H = H, H_tis = H_tis,
    cor = correlations,
    cor.p = ps,
    cross_cor = p_cross,
    iter = iter
  ))
}



run_BayesNMF.joint <- function(V0, V0_tis, W_init = NULL, sp = FALSE, 
                               alpha_proj, alpha_cross, n.iter = 2000, tol = 1e-7,
                               K = 10, alpha_tis = 0.1, n_reps = 10,
                               max_attempts = 30, random_seed = 1) {

  print(paste0("Running bNMF clustering procedure (", n_reps, " iterations)..."))
  print(sprintf("Using tolerance of %.2e!",tol))

  set.seed(random_seed)
  successful_runs <- list()
  attempt_counts <- 0

  while (length(successful_runs) < n_reps) {
    tryCatch({
      print(paste("Attempt:", length(successful_runs) + 1))
      res <- BayesNMF.L2EU.joint(
        V0 = V0, V0_tis = V0_tis, W_init = W_init, sp = sp, alpha_proj = alpha_proj,
        alpha_cross = alpha_cross, n.iter = n.iter, tol = tol, K = K, alpha_tis = alpha_tis
      )
      successful_runs[[length(successful_runs) + 1]] <- res
      attempt_counts <- 0  # Reset counter after a successful run
    }, error = function(e) {
      if (grepl("Maximum number of iterations reached without convergence.", e$message)) {
        stop("Maximum number of iterations reached without convergence. Stopping retries.")
      } else {
        # print(paste("Error occurred:", e$message))
        attempt_counts <- attempt_counts + 1  # Increment attempt counter for each retry
        if (attempt_counts >= max_attempts) {
          stop("Maximum number of attempts for a single run has been exceeded. Stopping.")
        }
      }
    })

    attempt_counts <- attempt_counts + 1  # Increment attempt counter for each retry

    if (attempt_counts >= max_attempts) {
      stop("Maximum number of attempts for a single run has been exceeded. Stopping.")
    }
  }
  return(successful_runs)
}

run_BayesNMF_with_fallback.joint <- function(V_traits, V_tissues,
                                       W_init,
                                       sp,
                                       alpha_proj,
                                       alpha_cross,
                                       tolerance,
                                       K,
                                       n.iter,
                                       alpha_tis,
                                       NbrReps,
                                       alpha_proj_step = 1,
                                       min_alpha_proj = 1,
                                       max_retries = 30,
                                       attempts = 30) {
  retries <- 0
  success <- FALSE
  result <- NULL

  while (!success && retries < max_retries && alpha_proj >= min_alpha_proj) {
    tryCatch({
      # Try running the function with the current alpha_proj
      result <- run_BayesNMF.joint(V0 = V_traits, V0_tis = V_tissues,
                             W_init = W_init,
                             sp = sp,
                             alpha_proj = alpha_proj,
                             alpha_cross = alpha_cross,
                             n.iter = n.iter,
                             tol = tolerance,
                             K = K,
                             alpha_tis = alpha_tis,
                             n_reps = NbrReps,
                             max_attempts = attempts)
      # If successful, set success to TRUE to exit the loop
      success <- TRUE
    }, error = function(e) {
      # If there's an error, reduce alpha_proj and retry
      cat("Error occurred with alpha_proj =", alpha_proj, "\n")
      cat("Error message:", e$message, "\n")
      alpha_proj <<- alpha_proj - alpha_proj_step # Update alpha_proj globally
      retries <<- retries + 1 # Update retries globally
      if (alpha_proj < min_alpha_proj) {
        cat("Minimum alpha_proj reached. Stopping retries.\n")
      } else {
        cat("Retrying with alpha_proj =", alpha_proj, "\n")
      }
    })
  }

  if (!success) {
    stop("Function failed after ", retries, " attempts due to non-convergence or failure to identify at least 2 clusters.")
  }

  return(result)
}

# ******************************************************
# EXTRA FUNCTIONS
# ******************************************************

# ******************************************************
# updates
# ******************************************************

update_H <- function(W, V, V.ap, H, lambda_proj, lambda_cross, alpha_proj, alpha_cross) {
  eps <- 1.e-50

  H_udpdate <- H * (t(W) %*% V) /
    (t(W) %*% V.ap +
    H * (alpha_proj * matrix(rep(1 / lambda_proj, ncol(V)), ncol = ncol(V)) +
    alpha_cross * matrix(rep(1 / lambda_cross, ncol(V)), ncol = ncol(V))) + eps)

  return(H_udpdate)
}


update_H_tis <- function(W, V_tis, V.ap_tis, H_tis, lambda_proj, lambda_cross, alpha_proj, alpha_cross) {
  eps <- 1.e-50

  H_tis_udpdate <- H_tis * (t(W) %*% V_tis) /
    (t(W) %*% V.ap_tis +
    H_tis * (alpha_proj * matrix(rep(1 / lambda_proj, ncol(V_tis)), ncol = ncol(V_tis)) +
    alpha_cross * matrix(rep(1 / lambda_cross, ncol(V_tis)), ncol = ncol(V_tis))) + eps)

  for (i in 1:ncol(H_tis_udpdate)) {
    max_val_index <- which.max(H_tis_udpdate[, i])
    H_tis_udpdate[-max_val_index, i] <- 0  # Set all other entries in the row to 0
  }

  return(H_tis_udpdate)
}


update_W_joint <- function(W, V, V.ap, V_tis, V.ap_tis, H, H_tis,
                           lambda_proj, lambda_cross, alpha_tis, alpha_proj, alpha_cross) {
  eps <- 1.e-50

  W * ((1 - alpha_tis) * V %*% t(H) + alpha_tis * V_tis %*% t(H_tis)) /
  ((1 - alpha_tis) * V.ap %*% t(H) + alpha_tis * V.ap_tis %*% t(H_tis) +
  W * (alpha_proj * t(matrix(rep( 1 / lambda_proj, nrow(V)), ncol = nrow(V))) +
  alpha_cross * t(matrix(rep(1 / lambda_cross, nrow(V)), ncol = nrow(V)))) + eps)
}

update_alpha_tis <- function(alpha_tis, V, V.ap, V_tis, V.ap_tis) {
  
  # Compute normalized matrices to compare them fairly
  V_norm <- V / sd(V); V.ap_norm <- V.ap / sd(V.ap)
  V_tis_norm <- V_tis / sd(V_tis); V.ap_tis_norm <- V.ap_tis / sd(V.ap_tis)

  # Correct reconstruction errors by their respective dimensionality
  err1 <- sum((V_norm - V.ap_norm)^2) / ncol(V)
  err2 <- sum((V_tis_norm - V.ap_tis_norm)^2) / ncol(V_tis)

  # Unconstrained update
  delta_alpha <- log(alpha_tis / (1 - alpha_tis)) + log(err2 / err1)
  
  # Apply sigmoid to constrain to [0, 1]
  alpha_tis_update <- exp(delta_alpha) / (1 + exp(delta_alpha))

  # alpha_tis_update <- alpha_tis * err2 / err1
  return(alpha_tis_update)
}


update_alpha_proj <- function(alpha_proj, V, H, V_tis, H_tis, V_ap, V_ap_tis, alpha_tis, max_alpha_proj) {
  # Compute projections
  V_proj <- V %*% t(H); V_proj_norm <- V_proj / sd(V_proj)
  V_tis_proj <- V_tis %*% t(H_tis); V_tis_proj_norm <- V_tis_proj / sd(V_tis_proj)
  
  # Compute alignment error
  alignment_error <- sum((V_proj_norm - V_tis_proj_norm)^2) / ncol(V_proj)
  
  # Compute reconstruction error
  reconstruction_error <- sum((1 - alpha_tis) * (V / sd(V) - V_ap / sd(V_ap))^2) / ncol(V) +
                          sum(alpha_tis * (V_tis / sd(V_tis) - V_ap_tis / sd(V_ap_tis))^2) / ncol(V_tis)
  
  # Update alpha_proj
  alpha_proj_update <- min(alpha_proj * (alignment_error / reconstruction_error), max_alpha_proj)
  
  return(alpha_proj_update)
}

# ******************************************************
# Simple gradient descent
# ******************************************************

# Function to obtain W from H_tissues and V_tissues
matrix_factorization_gradient_descent_nonneg_H <- function(V, W, sp = FALSE,
                                                           tol = 1e-5,
                                                           H_init = NULL,
                                                           max_iter = 10000,
                                                           seed = 1) {
  # plant seed
  set.seed(seed)

  T <- ncol(V)
  K <- ncol(W)

  # Initialize W with random values
  if (is.null(H_init)) {
    H <- matrix(runif(K * T), nrow = K, ncol = T)
  } else {
    H <- H_init
  }

  err <- 1
  iteration <- 1
  objective_value <- -99
  while (err > tol & iteration < max_iter) {

    objective_value_prev <- objective_value
    # Compute the difference between V and WH
    diff_V_WH <- V - W %*% H

    # Update W using the gradient
    #W <- W - learning_rate * gradient_W
    H <- H  * ((t(W) %*% V) / (t(W) %*% W %*% H + 1e-50))
    
    if(sp) {
      for (i in 1:ncol(H)) {
        max_val_index <- which.max(H[, i])
        H[-max_val_index, i] <- 0  # Set all other entries in the row to 0
      }
    }

    # Print the objective function value for monitoring
    objective_value <- sum(diff_V_WH^2)

    iteration <- iteration + 1
    err <- abs(objective_value - objective_value_prev)

   # cat(sprintf("Iteration %d: Objective Value = %.4f: err = %s\n", iteration,
   #             objective_value, err))
  }
  return(H)
}



matrix_factorization_gradient_descent_nonneg_W <- function(V, H, tol = 1e-5,
                                                           max_iter = 10000,
                                                           seed = 1) {
  # plant seed
  set.seed(seed)

  S <- nrow(V)
  K <- nrow(H)

  # Initialize W with random values
  W <- matrix(runif(S * K), nrow = S, ncol = K)

  err <- 1
  iteration <- 1
  objective_value <- -99
  while (err > tol & iteration < max_iter) {

    objective_value_prev <- objective_value
    # Compute the difference between V and WH
    diff_V_WH <- V - W %*% H

    # Update W using the gradient
    #W <- W - learning_rate * gradient_W
    V.ap <- W %*% H
    W <- W  * ((V %*% t(H)) / (V.ap %*% t(H) + 1e-50))

    # Print the objective function value for monitoring
    objective_value <- sum(diff_V_WH^2)

    iteration <- iteration + 1
    err <- abs(objective_value - objective_value_prev)

#    cat(sprintf("Iteration %d: Objective Value = %.4f: err = %s\n", iteration,
#                objective_value, err))
  }
  return(W)
}

# ******************************************************
# Merging clusters for Pleiotropic
# ******************************************************

merge_overlapping_clusters <- function(H_tissues, H_traits, W, threshold = 0.8) {
  # Compute similarity between clusters
  similarity <- as.matrix(proxy::simil(H_tissues, method = "cosine"))
  
  while (TRUE) {
    # Find the most similar pair of clusters
    max_sim <- max(similarity[lower.tri(similarity)])
    if (max_sim < threshold) break
    
    # Merge the two clusters
    merge_pair <- which(similarity == max_sim, arr.ind = TRUE)[1, ]
    new_cluster <- colMeans(H_tissues[merge_pair, ])
    new_cluster2 <- colMeans(H_traits[merge_pair, ])
    new_cluster3 <- rowMeans(W[, merge_pair])
    
    # Update the matrix
    H_tissues <- H_tissues[-merge_pair, , drop = FALSE]
    H_tissues <- rbind(H_tissues, new_cluster)

    H_traits <- H_traits[-merge_pair, , drop = FALSE]
    H_traits <- rbind(H_traits, new_cluster2)

    W <- W[, -merge_pair, drop = FALSE]
    W <- cbind(W, new_cluster3)
    # rownames(H_tissues)[nrow(H_tissues)] <- paste0("Merged_", merge_pair)
    
    # Recompute similarity
    similarity <- as.matrix(proxy::simil(H_tissues, method = "cosine"))
  }
  
  return(list(H_tissues = H_tissues, H_traits = H_traits, W = W))
}


# ******************************************************
# Recompute H_tissues for Pleiotropic to merge similar clusters
# ******************************************************

recompute_H_tis_pleio <- function(V_tis, H_tis, H, W, thres = 0.5) {
  H_tis <- matrix_factorization_gradient_descent_nonneg_H(V = V_tis, W = W, sp = FALSE)
  H_tis[H_tis < 1e-10] <- 0

  # Retain rows where at least one cluster is active
  valid_clusters <- rowSums(H_tis) >= 1e-10
  H <- H[valid_clusters, ]
  W <- W[, valid_clusters]
  H_tis <- H_tis[valid_clusters, ]

  # Merge similar clusters in terms of tissue profiles
  res <- merge_overlapping_clusters(H_tis, H, W, threshold = thres)
  H_tis <- res$H_tissues
  H <- res$H_traits
  W <- res$W

  # Enforce 1 cluster to tissue associations
  H_tis <- matrix_factorization_gradient_descent_nonneg_H(V = V_tis, W = W, H_init = H_tis, sp = TRUE)
  H_tis[H_tis < 1e-10] <- 0

  valid_clusters <- rowSums(H_tis) >= 1e-10
  H <- H[valid_clusters, ]
  W <- W[, valid_clusters]
  H_tis <- H_tis[valid_clusters, ]

  return(list(H_tis = H_tis, H = H, W = W))
}

# ******************************************************
# Rest
# ******************************************************

process_matrix_h_traits <- function(math, matw) {
  # Keep only rows where the sum is above or equal to 1.e-5
  math_upd <- math[rowSums(math) >= 1.e-10 & colSums(matw) >= 1.e-10, ]
  matw_upd <- matw[, rowSums(math) >= 1.e-10 & colSums(matw) >= 1.e-10]

  # Set all values in the matrix to be at least 1.e-10
  math_upd[math_upd < 1e-10] <- 0
  matw_upd[matw_upd < 1e-10] <- 0
  
  return(list(math = math_upd, matw = matw_upd))
}

process_matrix_h <- function(math, mathtis, matw) {
  # Keep only rows where the sum is above or equal to 1.e-5
  math_upd <- math[rowSums(math) >= 1.e-10 & rowSums(mathtis) >= 1.e-10 & colSums(matw) >= 1.e-10, ]
  mathtis_upd <- mathtis[rowSums(math) >= 1.e-10 & rowSums(mathtis) >= 1.e-10 & colSums(matw) >= 1.e-10, ]
  matw_upd <- matw[, rowSums(math) >= 1.e-10 & rowSums(mathtis) >= 1.e-10 & colSums(matw) >= 1.e-10]

  # Set all values in the matrix to be at least 1.e-10
  math_upd[math_upd < 1e-10] <- 0
  mathtis_upd[mathtis_upd < 1e-10] <- 0
  matw_upd[matw_upd < 1e-10] <- 0
  return(list(math = math_upd, mathtis = mathtis_upd, matw = matw_upd))
}

# function to make columns be the same between V and H
homogenize_columns <- function(V_traits, H_traits, V_tissues, H_tissues) {
  # List of matrices to process
  matrices <- list(V_traits, H_traits, V_tissues, H_tissues)
  # Corresponding column sets to intersect
  columns <- list(colnames(V_traits), colnames(H_traits),
                  colnames(V_tissues), colnames(H_tissues))

  # Intersect columns and update matrices
  for (i in seq(1, length(matrices), by = 2)) {
    common_cols <- intersect(columns[[i]], columns[[i + 1]])
    matrices[[i]] <- matrices[[i]][, common_cols]
    matrices[[i + 1]] <- matrices[[i + 1]][, common_cols]
  }

  # Return the updated matrices as a list
  return(list(V_traits = matrices[[1]],
              H_traits = matrices[[2]],
              V_tissues = matrices[[3]],
              H_tissues = matrices[[4]]))
}

# Function to load data from the given directory
load_data <- function(DIR, folder, Trait, glob_name, full_chr) {

  V_traits <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                              DIR, folder, Trait, "V_traits",
                              glob_name, full_chr))
  V_tissues <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                               DIR, folder, Trait, "V_tissues",
                               glob_name, full_chr))
  H_traits <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                              DIR, folder, Trait, "H_traits",
                              glob_name, full_chr))
  H_tissues <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                               DIR, folder, Trait, "H_tissues",
                               glob_name, full_chr))
  W_traits <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                              DIR, folder, Trait, "W_traits",
                              glob_name, full_chr))
  W_tissues <- readRDS(sprintf("%s/FIGURES%s/%s/%s_list_%s_loco_%s.rds",
                               DIR, folder, Trait, "W_tissues",
                               glob_name, full_chr))
  return(list(V_traits = V_traits, V_tissues = V_tissues,
              H_traits = H_traits, H_tissues = H_tissues,
              W_traits = W_traits, W_tissues = W_tissues))
}

calculate_covariances <- function(V_proj, V_tis_proj, H, H_tis) {
  myk <- which(rowSums(H) > 1e-10 & rowSums(H_tis) > 1e-10)
  covariances <- sapply(myk, function(j) {
    cov(V_proj[, j], V_tis_proj[, j])
  })
  return(covariances)
}

calculate_correlations <- function(V_proj, V_tis_proj, H, H_tis) {
  myk <- which(rowSums(H) > 1e-10 & rowSums(H_tis) > 1e-10)
  correlations <- sapply(myk, function(j) {
    cor(V_proj[, j], V_tis_proj[, j])
  })
  ps <- sapply(myk, function(j) {
    cor.test(V_proj[, j], V_tis_proj[, j])$p.val
  })
  return(list(correlations = correlations, p = ps))
}

calculate_cross_correlations <- function(V_proj, V_tis_proj, H, H_tis) {
  myk <- which(rowSums(H) > 1e-10 & rowSums(H_tis) > 1e-10)
  cross_correlations <- matrix(1, nrow = max(myk), ncol = max(myk))
  
  for (j in myk) {
    for (i in myk) {
      if (i != j) {
        # Compute correlation between V_proj[:, j] and V_tis_proj[:, i] and vice versa
        cross_correlations[j, i] <- max(cor.test(V_proj[, j], V_tis_proj[, i])$estimate, cor.test(V_tis_proj[, j], V_proj[, i])$estimate)
      }
    }
  }
  
  return(cross_correlations[myk, myk])
}

calculate_cross_correlationsp <- function(V_proj, V_tis_proj, H, H_tis) {
  myk <- which(rowSums(H) > 1e-10 & rowSums(H_tis) > 1e-10)
  cross_correlations <- matrix(1, nrow = max(myk), ncol = max(myk))
  
  for (j in myk) {
    for (i in myk) {
      if (i != j) {
        # Compute correlation between V_proj[:, j] and V_tis_proj[:, i] and vice versa
        cross_correlations[j, i] <- min(cor.test(V_proj[, j], V_tis_proj[, i])$p.value, cor.test(V_tis_proj[, j], V_proj[, i])$p.value)
      }
    }
  }
  
  return(cross_correlations[myk, myk])
}



# ******************************************************
# Functions for running, processing and saving BNMF models
# ******************************************************

# Define function to run the Bayesian NMF model
run_single_bNMF_model <- function(model_name, V_traits, V_tissues, Tolerance, k, NbrReps,
                           Alpha_proj = NULL, Alpha_cross = NULL, w_tis = NULL, attempts = 30) {
  
  # Initialize variables
  H_traits <- NULL
  try <- 1
  
  if (model_name %in% c("Pleiotropic", "Epigenomic")) {
    while ((is.null(dim(H_traits)) || dim(H_traits)[1] < 2) & try <= 10) {
      
      # Choose the appropriate model function
      if (model_name == "Pleiotropic") {
        model_func <- run_BayesNMF.trait
        input_args <- list(V0 = V_traits, W_init = NULL, sp = TRUE, tol = Tolerance, K = k,
                           n_reps = NbrReps, random_seed = try)
      } else {  # Epigenomic model
        model_func <- run_BayesNMF.tissue
        input_args <- list(V0 = V_tissues, W_init = NULL, sp = TRUE, tol = Tolerance, K = k,
                           n_reps = NbrReps, random_seed = try)
      }
      
      # Run Bayesian NMF for the current attempt
      bNMF_results <- do.call(model_func, input_args)
      res <- bNMF_results[[1]]
      
      # Process matrices based on model type
      if (model_name == "Pleiotropic") {
        W_traits <- process_matrix_h_traits(res$H, res$W)$matw
        H_traits <- process_matrix_h_traits(res$H, res$W)$math

        # Apply additional H_traits computation step for Epigenomic model
        if (!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2) { 
          H_tissues <- matrix_factorization_gradient_descent_nonneg_H(V = V_tissues + 1e-50, W = W_traits, sp = TRUE)
          H_traits <- H_traits[rowSums(H_tissues) >= 1.e-10, ]
        }
      } else {  # Epigenomic model
        W_traits <- process_matrix_h_traits(res$H, res$W)$matw
        H_tissues <- process_matrix_h_traits(res$H, res$W)$math
        
        # Apply additional H_traits computation step for Epigenomic model
        if (!is.null(dim(H_tissues)) && dim(H_tissues)[1] >= 2) { 
          H_traits <- matrix_factorization_gradient_descent_nonneg_H(V = V_traits + 1e-50, W = W_traits, sp = FALSE)
          H_tissues <- H_tissues[rowSums(H_traits) >= 1.e-10, ]
          H_traits <- H_traits[rowSums(H_traits) >= 1.e-10, ]
        }
      }
      
      try <- try + 1
    }
    
  } else if (model_name == "JPEP") {
    # JPEP does not require the while loop and additional processing steps
    model_func <- run_BayesNMF_with_fallback.joint
    input_args <- list(V_traits = V_traits, V_tissues = V_tissues, W_init = NULL, sp = TRUE, tolerance = Tolerance,
                       K = k, n.iter = 5000, alpha_proj = Alpha_proj, alpha_cross = Alpha_cross,
                       alpha_tis = w_tis, NbrReps = NbrReps, attempts = attempts)
    
    bNMF_results <- do.call(model_func, input_args)
  } else {
    stop("Invalid model name!")
  }
  
  return(bNMF_results)
}



process_results <- function(model_name, bNMF_results, V_traits, V_tissues) {
  
  # Initialize storage vectors
  cors_vec <- rep(NA, length(bNMF_results))
  ks <- rep(NA, length(bNMF_results))
  
  for (i in seq_len(length(bNMF_results))) {
    res <- bNMF_results[[i]]

    # Select the appropriate processing function
    if (model_name == "JPEP") {
      if (!(!is.null(dim(res$H)) && dim(res$H)[1] >= 2)) {
        next # Move to the next iteration
      }
      processed <- process_matrix_h(res$H, res$H_tis, res$W)
      H_tissues <- processed$mathtis
      H_traits <- processed$math
      if (!(!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2)) {
        next # Move to the next iteration
      }
    } else if (model_name == "Epigenomic") {
      if (!(!is.null(dim(res$H)) && dim(res$H)[1] >= 2)) {
        next # Move to the next iteration
      }
      processed <- process_matrix_h_traits(res$H, res$W)
      H_tissues <- processed$math
      # Skip iteration if fewer than 2 clusters are found
      if (!(!is.null(dim(H_tissues)) && dim(H_tissues)[1] >= 2)) {
        next # Move to the next iteration
      }
      H_traits <- matrix_factorization_gradient_descent_nonneg_H(V = V_traits, W = processed$matw, sp = TRUE)
      # Apply thresholding
      H_traits[H_traits < 1e-10] <- 0
      H_tissues <- H_tissues[rowSums(H_traits) >= 1.e-10, , drop = FALSE]
      processed$matw <- processed$matw[, rowSums(H_traits) >= 1.e-10, drop = FALSE]
      H_traits <- H_traits[rowSums(H_traits) >= 1.e-10, , drop = FALSE]
      if (!(!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2)) {
        next # Move to the next iteration
      }      
    } else if (model_name == "Pleiotropic") {
      if (!(!is.null(dim(res$H)) && dim(res$H)[1] >= 2)) {
        next # Move to the next iteration
      }
      processed <- process_matrix_h_traits(res$H, res$W)
      H_traits <- processed$math
      # Skip iteration if fewer than 2 clusters are found
      if (!(!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2)) {
        next # Move to the next iteration
      }
      H_tissues <- matrix_factorization_gradient_descent_nonneg_H(V = V_tissues, W = processed$matw, sp = TRUE)
      # Apply thresholding
      H_tissues[H_tissues < 1e-10] <- 0
      H_traits <- H_traits[rowSums(H_tissues) >= 1.e-10, , drop = FALSE]
      processed$matw <- processed$matw[, rowSums(H_tissues) >= 1.e-10, drop = FALSE]
      H_tissues <- H_tissues[rowSums(H_tissues) >= 1.e-10, , drop = FALSE]
      if (!(!is.null(dim(H_tissues)) && dim(H_tissues)[1] >= 2)) {
        next # Move to the next iteration
      }
    } else {
      stop("Invalid model name!")
    }

    if (!is.null(dim(processed$math)) && dim(processed$math)[1] >= 2) {
      
      # Compute correlation
      homogenized <- homogenize_columns(V_traits, H_traits, V_tissues, H_tissues)
      V_proj <- homogenized$V_traits %*% t(homogenized$H_traits)
      V_tis_proj <- homogenized$V_tissues %*% t(homogenized$H_tissues)
      cors <- calculate_correlations(V_proj, V_tis_proj, homogenized$H_traits, homogenized$H_tissues)
      
      cors_vec[i] <- mean(cors$correlations)
      ks[i] <- nrow(homogenized$H_traits)
    }
  }

  # Identify best clustering result
  if (!all(is.na(cors_vec))) {
    # Calculate the frequency of each cluster value k
    freq <- table(ks[!is.na(ks)])

    # Find the maximum frequency
    max_freq <- max(freq)

    # Get the numbers with the maximum frequency
    most_repeated <- as.numeric(names(freq[freq == max_freq]))

    # Choose the smallest number in case of a tie
    best_k <- min(most_repeated)

    best_it <- which(ks == best_k & cors_vec == max(cors_vec[which(ks == best_k)]))[1]

    res <- bNMF_results[[best_it]]

    # Apply the correct processing function again for the best result
    if (model_name == "JPEP") {
      processed <- process_matrix_h(res$H, res$H_tis, res$W)
      H_tissues <- processed$mathtis
      H_traits <- processed$math
      if (!(!is.null(dim(processed$math)) && dim(processed$math)[1] >= 2)) {
        stop("Less than 2 discovered clusters. No further analysis conducted.")
      }
    } else if (model_name == "Epigenomic") {
      processed <- process_matrix_h_traits(res$H, res$W)
      H_tissues <- processed$math
      if (!(!is.null(dim(H_tissues)) && dim(H_tissues)[1] >= 2)) {
        stop("Less than 2 discovered clusters. No further analysis conducted.")
      }
      H_traits <- matrix_factorization_gradient_descent_nonneg_H(V = V_traits, W = processed$matw, sp = TRUE)
      # Apply thresholding
      H_traits[H_traits < 1e-10] <- 0
      H_tissues <- H_tissues[rowSums(H_traits) >= 1.e-10, , drop = FALSE]
      processed$matw <- processed$matw[, rowSums(H_traits) >= 1.e-10, drop = FALSE]
      H_traits <- H_traits[rowSums(H_traits) >= 1.e-10, , drop = FALSE]
      if (!(!is.null(dim(H_tissues)) && dim(H_tissues)[1] >= 2)) {
        stop("Less than 2 discovered clusters. No further analysis conducted.")
      }
    } else if (model_name == "Pleiotropic") {
      processed <- process_matrix_h_traits(res$H, res$W)
      H_traits <- processed$math
      if (!(!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2)) {
        stop("Less than 2 discovered clusters. No further analysis conducted.")
      }
      H_tissues <- matrix_factorization_gradient_descent_nonneg_H(V = V_tissues, W = processed$matw, sp = TRUE)
      # Apply thresholding
      H_tissues[H_tissues < 1e-10] <- 0
      H_traits <- H_traits[rowSums(H_tissues) >= 1.e-10, , drop = FALSE]
      processed$matw <- processed$matw[, rowSums(H_tissues) >= 1.e-10, drop = FALSE]
      H_tissues <- H_tissues[rowSums(H_tissues) >= 1.e-10, , drop = FALSE]
      if (!(!is.null(dim(H_traits)) && dim(H_traits)[1] >= 2)) {
        stop("Less than 2 discovered clusters. No further analysis conducted.")
      }
    } else {
      stop("Invalid model name!")
    }

    save_list <- function(obj, name) {
      saveRDS(obj, sprintf("%s/%s_list_%s_loco_%s.rds", output_dir, name, glob_name, loco))
    }

    # output results
    return(list(V_trait = V_traits, V_tissue = V_tissues,
                W = processed$matw, H_trait = H_traits, H_tissue = H_tissues))
  } else {
    stop("Less than 2 discovered clusters. No further analysis conducted.")
  }
}





prepare_V_matrices <- function(Trait1, thres_pip, Thres_cor, Category, Add_annot,
                               Grouping, dir_out, Singlecell, loco, min_sec_pip, thres_aux) {

  # Retrieve epigenomic and pleiotropic data
  if (! Trait1 == "All_Metal_LDSC-CORR_Neff_comparison") {
    EPI <- getEpiPairAll(Trait1, pip_thres = thres_pip, thres_cor = Thres_cor, category = Category)
    epi_pair_all <- EPI$epi
    associated_traits <- EPI$assoc
    associated_traits_pip2 <- paste0(associated_traits, "_pip2")
    # Filter auxiliary traits based on PIP threshold
    valid_traits <- epi_pair_all[, lapply(.SD, function(i) sum(abs(i) > min_sec_pip)),
                                 .SDcols = associated_traits_pip2] >= thres_aux
    associated_traits <- associated_traits[valid_traits]
  } else {
    DIR <- "/n/groups/price/gaspard/PLEIOTROPY/EPIMAP/NOLAN/UKKBB_AND_NON_UKBB"
    epi_pair_all <- fread(sprintf("%s/MACS2/MULTI_TRAITS/%s_comparison.nosusie.Epimap.intersect.short.thres0.01.txt.gz", DIR, "All_Metal_LDSC-CORR_Neff"))
    associated_traits <- colnames(epi_pair_all)[10:ncol(epi_pair_all)]
    associated_traits <- associated_traits[!grepl("_pip2", associated_traits)]
    associated_traits_pip2 <- paste0(associated_traits, "_pip2")
    # Filter auxiliary traits based on PIP threshold
    valid_traits <- epi_pair_all[, lapply(.SD, function(i) sum(abs(i) > min_sec_pip)),
                                 .SDcols = associated_traits_pip2] >= thres_aux
    associated_traits <- associated_traits[valid_traits]
  }

  # Construct V_traits matrix
  V_traits <- as.matrix(epi_pair_all[, ..associated_traits])
  if (! Trait1 == "All_Metal_LDSC-CORR_Neff_comparison") {
    rownames(V_traits) <- paste(epi_pair_all[, CHR], epi_pair_all[, START], epi_pair_all[, STOP], sep = ";")
  } else {
    rownames(V_traits) <- epi_pair_all$RSID
    epi_pair_all[, CHR := paste0("chr", CHR)]
  } 
  # Retrieve SNP-tissue probabilistic assignments
  if (!Singlecell) {
    V_tissues <- get_Vtissues(Trait1, pip_thres = thres_pip, add_annot = Add_annot,
                              grouping = Grouping, tol = 1e-8, maxiter = 200, bckgd_correction = 1)
  } else {
  #  dir_vtis <- "/n/groups/price/gaspard/PLEIOTROPY/EPIMAP/NOLAN/V_tissues/singlecell"
  #  V_tissues <- readRDS(sprintf("%s/%s/V_tissues_tol1e-10_groupingTRUE", dir_vtis, Trait1))
     if (! Trait1 == "All_Metal_LDSC-CORR_Neff_comparison") {
       V_tissues <- get_Vtissues(Trait1, pip_thres = thres_pip, add_annot = Add_annot, singlecell = TRUE,
                                 grouping = FALSE, tol = 1e-8, maxiter = 200, bckgd_correction = 1)
     } else {
       V_tissues <- get_Vtissues("All_Metal_LDSC-CORR_Neff", pip_thres = thres_pip, add_annot = Add_annot, singlecell = TRUE,
                                 grouping = FALSE, tol = 1e-8, maxiter = 200, bckgd_correction = 1)
     }
  }

  # Restrict V_tissues to S-LDSC derived causal tissues if applicable
  if (grepl("LDSC-prior", dir_out)) {
    cor <- data.table(aux_trait = associated_traits)
    causal_tis_dt <- get_causal_tissues(Trait1, cor, thres_tissue_sig = 0.2, colnames(V_tissues))
    sig_tissues_pt <- causal_tis_dt$tissues
    if (length(sig_tissues_pt) <= 1) sig_tissues_pt <- NULL

    V_tissues <- get_Vtissues(Trait1, pip_thres = thres_pip, add_annot = Add_annot, grouping = Grouping, sig_tis = sig_tissues_pt, tol = 1e-8, maxiter = 200, bckgd_correction = 1)
  }

  # LOCO: Remove specified chromosome
  if (loco %in% 1:22 && epi_pair_all[CHR == paste0("chr", loco), .N] > 1) {
    epi_pair_all <- epi_pair_all[!(CHR == paste0("chr", loco)), ]
    V_traits <- as.matrix(epi_pair_all[, ..associated_traits])
    rownames(V_traits) <- paste(epi_pair_all[, CHR], epi_pair_all[, START], epi_pair_all[, STOP], sep = ";")
    common_snps <- intersect(rownames(V_traits), rownames(V_tissues))
    V_traits <- V_traits[common_snps, , drop = FALSE]
    V_tissues <- V_tissues[common_snps, , drop = FALSE]    
    # V_tissues <- V_tissues[which(rownames(V_tissues) %in% rownames(V_traits)), ]
    V_tissues <- V_tissues[rownames(V_traits), ]
  } else if (loco %in% 1:22) {
    stop("Requested LOCO approach but less than 2 SNPs on chromosome.")
  }

  # Ensure consistent ordering across matrices
  common_snps <- intersect(rownames(V_traits), rownames(V_tissues))
  V_traits <- V_traits[common_snps, , drop = FALSE]
  V_tissues <- V_tissues[common_snps, , drop = FALSE]
  V_tissues <- V_tissues[rownames(V_traits), ]

  # Extend auxiliary traits to include positive and negative pleiotropic components
  V_traits <- cbind(
    pmax(V_traits, 0),
    pmax(-V_traits, 0)
  )
  colnames(V_traits) <- c(paste0(colnames(epi_pair_all[, ..associated_traits]), "_pos"),
                          paste0(colnames(epi_pair_all[, ..associated_traits]), "_neg"))

  # Filter out rows where all values are zero
  non_zero_rows <- rowSums(V_traits) > 0
  V_traits <- V_traits[non_zero_rows, , drop = FALSE]
  V_tissues <- V_tissues[non_zero_rows, , drop = FALSE]

  # Extract and align PIPs
  epi_pair_all <- epi_pair_all[paste(epi_pair_all$CHR, epi_pair_all$START, epi_pair_all$STOP, sep = ";") %in% rownames(V_traits), ]
  
  # Reorder epi_pair_all to match the order of V_traits
  epi_pair_all <- epi_pair_all[match(rownames(V_traits), paste(epi_pair_all$CHR, epi_pair_all$START, epi_pair_all$STOP, sep = ";")), ]

  pips <- epi_pair_all[, PIP]

  return(list(V_traits = V_traits, V_tissues = V_tissues, pips = pips))
}


# load Pleiotropic and Epigenomic outputs to define the maximum number of clusters to identify
cluster_prior <- function(Model, base_name, glob_name, Trait1) {

  epi_name <- "Epigenomic_sparse"
  pleio_name <- "Udler_sparse"

  DIR_EPI <- sprintf("%s%s", base_name, epi_name)
  DIR_PLEIO <- sprintf("%s%s", base_name, pleio_name)
    
  data_out <- load_data(DIR, DIR_EPI, Trait1, glob_name, "ALL")
  W_traits_epi <- data_out$W_traits[[1]]
  n1 <- ncol(W_traits_epi)

  data_out <- load_data(DIR, DIR_PLEIO, Trait1, glob_name, "ALL")
  W_traits_pleio <- data_out$W_traits[[1]]
  n2_0 <- ncol(W_traits_pleio)
  H_traits_pleio <- data_out$H_traits[[1]]
  H_tissues_pleio <- data_out$H_tissues[[1]]
  V_traits_pleio <- data_out$V_traits[[1]]
  V_tissues_pleio <- data_out$V_tissues[[1]]

  # Recompute H_tissues with sparsity condition (specific to Pleiotropic analysis)
  new_mat <- recompute_H_tis_pleio(V_tissues_pleio, H_tissues_pleio,
                                   H_traits_pleio, W_traits_pleio)
  H_tissues_pleio <- new_mat$H_tis
  H_traits_pleio <- new_mat$H
  W_traits_pleio <- new_mat$W

  n2 <- ncol(W_traits_pleio)
  if (Model == "JPEP") {
    return(max(c(n1, n2)))
  } else if (Model == "Pleiotropic") {
    return(n1)
  } else {
    return(n2_0)
  }
}

