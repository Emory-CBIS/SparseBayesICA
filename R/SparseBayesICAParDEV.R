#---- Main Function
# Y <- Q x V x N
# initial_values should contain A, S0, and beta








#---- Main Function
# Y <- Q x V x N
# initial_values should contain A, S0, and beta
SparseBayesICAParDEV <- function(cl, Y, X, initial_values,
                              nburnin = 10, nmcmc = 10,
                              concentration = 1.0,
                              print_freq = 100){


  # Initial Setup
  Q <- dim(initial_values$S0)[2]
  V <- dim(Y)[1]
  N <- dim(X)[1]
  P <- dim(X)[2]

  # DPM hyperparameters
  prior_shape <- 1.0
  prior_scale <- 1.0

  # Angles for the mixing matrix update
  angle_range = seq(from = 0.0, to = 2*pi, length.out = 100)
  sin_theta <- sin(angle_range)
  cos_theta <- cos(angle_range)

  # Storage for suff stat for mixing matrix update
  YSiall      <- array(0, dim = c(Q*N, Q))
  YSiOmega    <- array(0, dim = c(Q*N, Q))
  Aall        <- initial_values$A + 0

  S0_DPM = initial_values$S0_DPM
  cluster_mappings = rep(0, length(S0_DPM$n_in_cluster))
  n_in_cluster <- S0_DPM$n_in_cluster

  sigma_sq_q <- initial_values$sigma_sq_q + 0


  tau_sum = matrix(0, nrow = Q, ncol = P)
  tau_sq = initial_values$Horseshoe$tau_sq + 0
  xi     = initial_values$Horseshoe$xi + 0

  if (min(S0_DPM$cluster_membership) == 1){
    S0_DPM$cluster_membership = S0_DPM$cluster_membership - as.integer(1)
  }

  # https://stackoverflow.com/questions/22739876/how-to-export-objects-to-parallel-clusters-within-a-function-in-r
  export_to_cluster(cl, 'nmcmc', nmcmc)

  # Storage - posterior quantities
  Ai_posterior_mean         <- array(0, dim = c(Q*N, Q))
  sigma_sq_q_posterior_mean <- rep(0, Q)
  tau_sq_posterior_mean     <- matrix(0, nrow = Q, ncol = P)

  # Verify input dimensions all match

  # Initialize object to keep track of computation time
  timing <- create_sparsebayes_timing_log()

  # Obtain the cluster splits
  n_node <- length(cl)
  voxel_combinations   <- partition_into_chunks(n_node, V)
  subject_combinations <- partition_into_chunks(n_node, N)
  subject_combinations_NQ <- list()
  for (ic in 1:n_node){
    inds <- ((subject_combinations[[ic]][1]-1)*Q+1):(tail(subject_combinations[[ic]], 1) * Q)
    subject_combinations_NQ[[ic]] <- inds
  }

  # Initialize nodes
  initialize_sparsebayes_nodes(cl,
                               voxel_combinations,
                               subject_combinations_NQ,
                               Y, X,
                               initial_values,
                               angle_range, sin_theta, cos_theta,
                               S0_DPM)

  actual_cluster_count <- 5

  overall_start_time <- Sys.time()
  for (imcmc in 1:(nburnin + nmcmc)){

    print("")
    get_perf_temp(Y,
                  Aall,
                  clusterEvalQ(cl, S0)[[1]],
                  clusterEvalQ(cl, Beta)[[1]],
                  X)

    if (round(imcmc / print_freq) == (imcmc / print_freq) ){
      print(paste0("ITER ", imcmc))
      cat("Means are:")
      print(S0_DPM$miu_h[1:actual_cluster_count])
      cat("Variances are:")
      print(S0_DPM$sigma_sq_h[1:actual_cluster_count])
      cat("Cluster counts are:")
      print(n_in_cluster[1:actual_cluster_count])
    }

    # Store prior precision
    Omega <- clusterEvalQ(cl, Omega <- diag(1.0 / sigma_sq_q))[[1]]

    #
    # Mixing Matrix Update
    #

    # Monitor Mixing Matrix Update Time
    mixing_matrix_start_time = Sys.time()

    # Get the required quantities for the mbvmf update
    clusterEvalQ(cl, calculate_YiSit_inplace(YSi, QxQStore, Y, Si, Q))
    elapsed_time <- Sys.time() - mixing_matrix_start_time
    timing[which(timing$Name == "YiSi Update"),]$Value =
      timing[which(timing$Name == "YiSi Update"),]$Value + as.numeric(elapsed_time)

    start_time <- Sys.time()
    YSiall[,]      <- Reduce(`+`, clusterEvalQ(cl, YSi))
    matmul_inplace(YSiOmega, YSiall, Omega)
    elapsed_time   <- Sys.time() - start_time
    timing[which(timing$Name == "YSi Reduction"),]$Value =
      timing[which(timing$Name == "YSi Reduction"),]$Value + as.numeric(elapsed_time)


    if (any(is.na(YSiOmega)) || any(is.na(YSiall))){
      stop("Missing values in Y x Si")
    }

    start_time <- Sys.time()
    # Export each piece to the corresponding node
    negYSiSigInvover2 <- -Omega/2
    update_cluster_matrix_variable(cl, 'negYSiSigInvover2', negYSiSigInvover2)
    update_cluster_matrix_variable_rowspec(cl, 'YSiOmegaChunk', YSiOmega, 'mixing_matrix_node_row_inds')

    elapsed_time <- Sys.time() - start_time
    timing[which(timing$Name == "mixing matrix partition export"),]$Value =
      timing[which(timing$Name == "mixing matrix partition export"),]$Value + as.numeric(elapsed_time)


    start_time <- Sys.time()
    clusterEvalQ(cl, gibbs_sample_mixing_matrix_inplace(A,
                                               YSiOmegaChunk,
                                               negYSiSigInvover2,
                                               YYt,
                                               angle_range,
                                               sin_theta,
                                               cos_theta,
                                               column_reduced_mixing_matrix,
                                               null_space))
    elapsed_time <- Sys.time() - start_time
    timing[which(timing$Name == "mixing matrix sampler"),]$Value =
      timing[which(timing$Name == "mixing matrix sampler"),]$Value + as.numeric(elapsed_time)

    # DELETE ME
    #Aallprev <- Aall + 0

    start_time <- Sys.time()
    # Collect updated samples from all nodes
    Aall[,] = abind::abind(clusterEvalQ(cl, A), along=1)
    elapsed_time <- Sys.time() - start_time
    timing[which(timing$Name == "mixing matrix collection"),]$Value =
      timing[which(timing$Name == "mixing matrix collection"),]$Value + as.numeric(elapsed_time)

    print("")
    print("After update A")
    get_perf_temp(Y,
                  Aall,
                  clusterEvalQ(cl, S0)[[1]],
                  clusterEvalQ(cl, Beta)[[1]],
                  X)

#
#     # DELETE ME
#     # this runs a check for any label switching
#     label_swaps <- 0
#     for (i in 1:N){
#       iS <- (i-1)*Q + 1
#       iE <- i*Q
#       column_correlation <- compare_2d_map_correlations(Aall[iS:iE,], Aallprev[iS:iE,])
#       # Top correlation should also be diagonal, otherwise register a label swap
#       if (any(column_correlation$sorted_order != 1:Q)){
#         label_swaps <- label_swaps + 1
#       }
#     }
#
#     if (label_swaps > 0){
#       print(paste0("LABELS WERE SWITCHED IN ", label_swaps, " INSTANCES"))
#       stop("lsw")
#     }

    # Sync copies of the mixing matrix across all nodes
    # The version with all subjects is stored as Ai, not A
    start_time <-Sys.time()
    #export_to_cluster(cl, 'Ai', Aall)
    update_cluster_matrix_variable(cl, 'Ai', Aall)
    elapsed_time <- Sys.time() - start_time
    timing[which(timing$Name == "mixing matrix export"),]$Value =
      timing[which(timing$Name == "mixing matrix export"),]$Value + as.numeric(elapsed_time)

    # Finish Monitoring Mixing Matrix Update Time
    elapsed_time <- Sys.time() - mixing_matrix_start_time
    timing[which(timing$Name == "Overall Mixing Matrix Time"),]$Value =
      timing[which(timing$Name == "Overall Mixing Matrix Time"),]$Value + as.numeric(elapsed_time)

    #
    # AtY Bookkeeping - is Ai' Yi
    #

    AtY_start_time <- Sys.time()

    clusterEvalQ(cl, calculate_AitYi_inplace(AtY, VxQStore, Ai, Y, Q))

    elapsed_time <- Sys.time() - AtY_start_time
    timing[which(timing$Name == "At x Y Time"),]$Value =
      timing[which(timing$Name == "At x Y Time"),]$Value + as.numeric(elapsed_time)
    timing[which(timing$Name == "Overall Variable Management Time"),]$Value =
      timing[which(timing$Name == "Overall Variable Management Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample DPM Quantities
    #

    cluster_membership_start_time <- Sys.time()

    # Get the minimum value of u (across all nodes)
    min_u = min(unlist(clusterEvalQ(cl, min(u))))

    # Sample stick breaking weights
    H <- which(n_in_cluster == 0)[1] - 1
    if (is.na(H)){H = length(n_in_cluster)}
    n_cluster      <- H + 10
    if (n_cluster > length(n_in_cluster)){
      print("ADJUSTING STORAGE SIZE!")
      nAdd <- n_cluster - length(n_in_cluster)
      n_in_cluster <- c(n_in_cluster, rep(0, nAdd))
      S0_DPM$miu_h      <- c(S0_DPM$miu_h,  rnorm(nAdd))
      S0_DPM$sigma_sq_h <- c(S0_DPM$sigma_sq_h,  invgamma::rinvgamma(nAdd, 1, 1))
      cluster_mappings <- c(cluster_mappings, rep(0, nAdd))
    }
    n_remaining    <- Q*V - cumsum(n_in_cluster[1:n_cluster])
    draws          <- rbeta(n_cluster, 1 + n_in_cluster[1:n_cluster], concentration + n_remaining)
    stick_remaining <- cumprod(c(1, 1 - draws))
    S0_DPM$miu_h[(H+1):n_cluster]      <- rnorm(10)
    S0_DPM$sigma_sq_h[(H+1):n_cluster] <- invgamma::rinvgamma(10, 1, 1)

    export_to_cluster(cl, 'mu_h', S0_DPM$miu_h)
    export_to_cluster(cl, 'sigma_sq_h', S0_DPM$sigma_sq_h)
    S0_DPM$sticks[1:n_cluster] <- draws * stick_remaining[-length(draws)]
    export_to_cluster(cl, 'n_in_cluster', n_in_cluster)

    # Draw any additional required clusters
    if (min_u < S0_DPM$sticks[n_cluster]){
      print("NEED TO DRAW MORE STICKS!")
    }

    # Update sticks on clusters
    #export_to_cluster(cl, 'sticks', S0_DPM$sticks)
    #stick_temp <- S0_DPM$sticks
    update_cluster_vector_variable(cl, 'sticks', S0_DPM$sticks)


    # Sample U
    clusterEvalQ(cl, sample_u(u, sticks, cluster_memberships))

    # Sample cluster memberships
    export_to_cluster(cl, 'H', n_cluster)
    #export_to_cluster(cl, 'n_in_cluster', n_in_cluster)
    update_cluster_vector_variable(cl, 'n_in_cluster', n_in_cluster)

    clusterEvalQ(cl, sample_cluster_membership_indicators(cluster_memberships,
                                                          n_in_cluster,
                                                          S0, u, sticks, mu_h,
                                                          sigma_sq_h, sigma_sq_q, H))


    # Summarize cluster counts across nodes
    n_in_cluster <- Reduce(`+`, clusterEvalQ(cl, n_in_cluster))
    update_cluster_vector_variable(cl, 'n_in_cluster', n_in_cluster)

    # Cleanup
    generate_cluster_mapping_vector_inplace(cluster_mappings, n_in_cluster)
    export_to_cluster(cl, 'mapping', cluster_mappings)
    #update_cluster_vector_variable(cl, 'mapping', cluster_mappings)

    actual_cluster_count <- sum(cluster_mappings > 0)
    if (cluster_mappings[1] == 0){actual_cluster_count <- actual_cluster_count + 1}

    export_to_cluster(cl, 'actual_cluster_count', actual_cluster_count)

    clusterEvalQ(cl, cleanup_cluster_ordering_inplace(cluster_memberships, mu_h,
                                                      sigma_sq_h, n_in_cluster,
                                                      mapping, actual_cluster_count))

    # Reorder the mu, sigma for the DPM
    S0_DPM$miu_h <- clusterEvalQ(cl[1], mu_h)[[1]]
    S0_DPM$sigma_sq_h <- clusterEvalQ(cl[1], sigma_sq_h)[[1]]
    n_in_cluster <- clusterEvalQ(cl[1], n_in_cluster)[[1]]

    elapsed_time <- Sys.time() - cluster_membership_start_time
    timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value =
      timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample S0 and Beta
    #

    spatial_map_start_time <- Sys.time()

    export_to_cluster(cl, 'mu_h_div_sigma_sq_h', S0_DPM$miu_h / S0_DPM$sigma_sq_h)

    cat("\n")
    L1magBeta <- sum(abs((clusterEvalQ(cl, Beta)[[1]])))
    L1magS0 <- sum(abs((clusterEvalQ(cl, S0)[[1]])))
    print("Before update  spatial maps")
    print(paste0("Magnitude of S0: ", L1magS0, " and magnitude of Beta: ", L1magBeta))

    clusterEvalQ(cl, sample_spatial_maps_inplace(S0, Beta, posterior_mean,
                                                 posterior_variance,
                                                prior_precision, draw,
                                                AtY, X_with_int, XtX, sigma_sq_q,
                                                tau_sq, lambda_sq,
                                                mu_h_div_sigma_sq_h, sigma_sq_h,
                                                cluster_memberships,
                                                Ip))


    if ( any(unlist(clusterEvalQ(cl, is.na(Beta)))) ){
     stop("NA beta was obtained")
    }


    cat("\n")
    L1magBeta <- sum(abs((clusterEvalQ(cl, Beta)[[1]])))
    L1magS0 <- sum(abs((clusterEvalQ(cl, S0)[[1]])))
    print("After update  spatial maps")
    print(paste0("Magnitude of S0: ", L1magS0, " and magnitude of Beta: ", L1magBeta))
    get_perf_temp(Y,
                  Aall,
                  clusterEvalQ(cl, S0)[[1]],
                  clusterEvalQ(cl, Beta)[[1]],
                  X)



    if ( any(unlist(clusterEvalQ(cl, is.na(S0)))) ){
      stop("NA S0 was obtained")
    }

    elapsed_time <- Sys.time() - spatial_map_start_time
    timing[which(timing$Name == "Overall Spatial Map Time"),]$Value =
      timing[which(timing$Name == "Overall Spatial Map Time"),]$Value + as.numeric(elapsed_time)


    #
    # Sample DPM Hyperparameters
    #

    DPM_hyperparameter_time <- Sys.time()

    # Evaluate the "sufficient statistics"
    clusterEvalQ(cl, calculate_DPM_hyperparameter_suffstats_inplace(sum_Ssq_qh,
                                                                   Sbar_qh,
                                                                   n_in_cluster_qh,
                                                                   S0,
                                                                   cluster_memberships))

    # Unite the sufficient statistics across the different nodes
    sum_Ssq_qh_total      <- Reduce(`+`, clusterEvalQ(cl, sum_Ssq_qh) )
    n_in_cluster_qh_total <- Reduce(`+`, clusterEvalQ(cl, n_in_cluster_qh) )

    # Account for chance that none in cluster -> can cause nan
    n_in_cluster_qh_total_temp <- n_in_cluster_qh_total
    n_in_cluster_qh_total_temp[n_in_cluster_qh_total_temp == 0] <- 1
    Sum_S_qh              <- Reduce(`+`, clusterEvalQ(cl, Sbar_qh) )
    Sbar_qh_total         <- Sum_S_qh / n_in_cluster_qh_total_temp # also not scaled by variance yet

    # Now obtain the S - Sbar squared term - specific to each node
    # TODO - get rid of this quantity
    DPM_hyper_sse_unscaled <- (Sum_S_qh^2 - 2*Sum_S_qh*Sbar_qh_total + Sbar_qh_total^2)

    sample_DPM_hyperparameters_inplace(S0_DPM$miu_h,
                                       S0_DPM$sigma_sq_h,
                                       sum_Ssq_qh_total,
                                       Sbar_qh_total,
                                       DPM_hyper_sse_unscaled,
                                       n_in_cluster_qh_total,
                                       n_in_cluster,
                                       sigma_sq_q,
                                       prior_shape,
                                       prior_scale,
                                       actual_cluster_count)

    #print("cheating")
    #S0_DPM$miu_h[1] <- 0


    # Push new cluster parameters to each worker process
    export_to_cluster(cl, 'mu_h', S0_DPM$miu_h)
    export_to_cluster(cl, 'sigma_sq_h', S0_DPM$sigma_sq_h)

    # # DELETE
    # # getting a calculation for how many have each cluster available:
    # uu <- clusterEvalQ(cl, u)[[1]]
    # nav <- rep(0, 20)
    # for (h in 1:20){
    #   nav[h] = sum(uu < S0_DPM$sticks[h])
    # }
    # mu_param <- S0_DPM$miu_h[1:20]
    # sig_param <- S0_DPM$sigma_sq_h[1:20]
    # mu_est <- rep(NA, 20)
    # sig_est <- rep(NA, 20)
    # myS0 <- clusterEvalQ(cl, S0)[[1]]
    # myCM <- clusterEvalQ(cl, cluster_memberships)[[1]]
    # for (h in 1:20){
    #   if (n_in_cluster[h] > 0){
    #     mu_est[h] <- mean(myS0[which(myCM == (h-1))])
    #     sig_est[h] <- var(myS0[which(myCM == (h-1))])
    #   }
    # }
    # print(cbind( S0_DPM$sticks[1:20], n_in_cluster[1:20], nav, mu_param, sig_param, mu_est, sig_est))


    if (any(is.na(S0_DPM$miu_h))){
      stop("NA cluster mean obtained")
    }

    # Update concentration parameter
    phi <- rbeta(1, concentration + 1, Q*V)
    odds <- (4 + H - 1) / (Q*V * (4.0 - log(phi)))

    if (runif(1, 0,1) < (odds/(odds+1)) ){
      concentration <- rgamma(1, 4.0 + H, scale = 1.0 / (4.0 - log(phi)))
    } else {
      concentration <- rgamma(1, 4.0 + H - 1.0, scale = 1.0 / (4.0 - log(phi)))
    }



    elapsed_time <- Sys.time() - DPM_hyperparameter_time
    timing[which(timing$Name == "Overall DPM Hyperparameter Time"),]$Value =
      timing[which(timing$Name == "Overall DPM Hyperparameter Time"),]$Value + as.numeric(elapsed_time)

    #
    # Update Errors predicted subject level maps
    #

    update_eij_sij_start_time <- Sys.time()

    # Update the errors
    clusterEvalQ(cl, evaluate_error_inplace(Ei, Si, AtY, S0, Beta, X_with_int))

    elapsed_time <- Sys.time() - update_eij_sij_start_time
    timing[which(timing$Name == "Eij and Sij Time"),]$Value =
      timing[which(timing$Name == "Eij and Sij Time"),]$Value + as.numeric(elapsed_time)
    timing[which(timing$Name == "Overall Variable Management Time"),]$Value =
      timing[which(timing$Name == "Overall Variable Management Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample Variance Terms
    #

    overall_variance_terms_time <- Sys.time()

    # Evaluate SSE
    clusterEvalQ(cl, evaluate_sse_inplace(sse, Ei))

    # Evaluate the regression scale term
    clusterEvalQ(cl, evaluate_reg_scale_inplace(reg_scale, S0, Beta, mu_h, sigma_sq_h,
                                               tau_sq, cluster_memberships, lambda_sq))

    # Aggregate across worker processes
    sse_total       <- Reduce(`+`, clusterEvalQ(cl, sse) )
    reg_scale_total <- Reduce(`+`, clusterEvalQ(cl, reg_scale) )

    # Draw new variance parameters
    shape      = (N * V + (P+1) * V) / 2.0
    scale      = sse_total / 2.0 + reg_scale_total / 2.0
    sigma_sq_q = invgamma::rinvgamma(Q, shape, rate = scale)
    #sigma_sq_q = sigma_sq_q / sigma_sq_q
    export_to_cluster(cl, 'sigma_sq_q', sigma_sq_q)

    # Check if failed
    if (any(is.na(sigma_sq_q))){
      print("Sigma squared h")
      print(S0_DPM$sigma_sq_h)
      print("n in cluster")
      print(n_in_cluster)
      print("SSE")
      print(sse_total)
      print("Reg Scale term:")
      print(reg_scale_total)
      stop("NA variance obtained")
    }

    # Local Shrinkage Parameters
    clusterEvalQ(cl, sample_local_shrinkage_parameters_inplace(lambda_sq,
                                                               nu,
                                                               cauchy_mixing,
                                                               cauchy_mixing_prior,
                                                               Beta,
                                                               tau_sq,
                                                               sigma_sq_q))

    if (any(unlist(clusterEvalQ(cl, sum(is.na(lambda_sq)))))){
      stop("NA Lambda obtained")
    }

    # Global (Group) Shrinkage Parameters
    clusterEvalQ(cl, calculate_tausum_inplace(tau_sum, Beta, lambda_sq, sigma_sq_q))

    tau_sum[,] = Reduce(`+`, clusterEvalQ(cl, tau_sum) )

    for (q in 1:Q){
      for (p in 1:P){
        tau_sq[q, p] = max(invgamma::rinvgamma(1,  (V+1)/2, rate = 1/xi[q,p] + tau_sum[q,p] ), 1e-20)
        xi[q, p]     =  invgamma::rinvgamma(1,  1, rate = 1 + 1/tau_sq[q,p] )
      }
    }

    # print("Tau")
    # print(tau_sq)

    if (any(is.na(tau_sq))){
      print("Tau sum")
      print(tau_sum)
      print("Tau sum on nodes")
      print(clusterEvalQ(cl, tau_sum))
      print("Lambda")
      print(clusterEvalQ(cl, lambda_sq))
      stop("NA Tau-squared obtained")
    }

    export_to_cluster(cl, 'tau_sq', tau_sq)

    elapsed_time <- Sys.time() - overall_variance_terms_time
    timing[which(timing$Name == "Overall Variance Term Time"),]$Value =
      timing[which(timing$Name == "Overall Variance Term Time"),]$Value + as.numeric(elapsed_time)


    #
    # Update Posterior Tracking
    #

    if (imcmc > nburnin){

      sample_storage_time <- Sys.time()

      Ai_posterior_mean[,] = Ai_posterior_mean +  1 / nmcmc * Aall

      sigma_sq_q_posterior_mean = sigma_sq_q_posterior_mean + 1/nmcmc * sigma_sq_q

      tau_sq_posterior_mean = tau_sq_posterior_mean + 1/nmcmc * tau_sq

      clusterEvalQ(cl,
                   update_spatial_map_posterior_tracking_inplace(S0_posterior_mean,
                                                                 Beta_posterior_mean,
                                                                 Beta_ngt0,
                                                                 lambda_sq_posterior_mean,
                                                                 S0, Beta, lambda_sq, nmcmc))



      elapsed_time <- Sys.time() - sample_storage_time
      timing[which(timing$Name == "Overall Sample Storage Time"),]$Value =
        timing[which(timing$Name == "Overall Sample Storage Time"),]$Value + as.numeric(elapsed_time)




    }

  }

  overall_run_time <- difftime(Sys.time(), overall_start_time, units = "secs")
  timing[which(timing$Name == "Overall MCMC Time"),]$Value = as.numeric(overall_run_time)

  # Reshape variables to more easily understood format
  Ai_3D <- array(0, dim = c(Q, Q, N))
  for (i in 1:N){
    Ai_temp  <- Ai_posterior_mean[ ((i-1)*Q+1):(i*Q), ]
    Ai_ortho <- Ai_temp %*% sqrtm(solve(t(Ai_temp)%*%Ai_temp));
    Ai_3D[,,i] <- Ai_ortho
  }

  S0_PM   <- abind::abind(clusterEvalQ(cl, S0_posterior_mean), along=1)
  Beta_PM <- abind::abind(clusterEvalQ(cl, Beta_posterior_mean), along=1)
  lambda_sq_PM <- abind::abind(clusterEvalQ(cl, lambda_sq_posterior_mean), along=1)

  # Get p-values based on credible intervals
  Beta_prop_gt0 <- abind::abind(clusterEvalQ(cl, Beta_ngt0), along=1) / nmcmc
  Beta_pvals <- 1 - Beta_prop_gt0
  Beta_pvals[Beta_pvals > 0.5] <- 1.0 - Beta_pvals[Beta_pvals > 0.5]
  Beta_pvals <- 2.0 * Beta_pvals # two-sided correction

  # Reorder Betas and P-values into 3d array, where 3rd dimension indexes covariates
  Beta_PM_reordered <- array(0, dim = c(V, Q, P))
  lambda_sq_PM_reordered <- array(0, dim = c(V, Q, P))
  Beta_pvalues_reordered <- array(0, dim = c(V, Q, P))
  for (q in 1:Q){
    index_set <- ((q-1)*P+1):(q*P)
    Beta_PM_reordered[,q,]      <- Beta_PM[, index_set]
    Beta_pvalues_reordered[,q,] <- Beta_pvals[, index_set]
    lambda_sq_PM_reordered[,q,] <- lambda_sq_PM[, index_set]
  }

  return(list(
    Ai_posterior_mean = Ai_3D,
    S0_posterior_mean = S0_PM,
    Beta_posterior_mean = Beta_PM_reordered,
    Beta_pvals = Beta_pvalues_reordered,
    sigma_sq_q_posterior_mean = sigma_sq_q_posterior_mean,
    tau_sq_posterior_mean = tau_sq_posterior_mean,
    lambda_sq_posterior_mean = lambda_sq_PM_reordered,
    timing = timing)
  )

}

