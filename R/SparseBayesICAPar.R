#' Internal function running the MCMC when a parallel environment is used (cl)
#' this is called by fit_SparseBayesICA, and is not a public function
#' @keywords internal
SparseBayesICAPar <- function(cl, Y, X, initial_values,
                              nburnin = 5000, nmcmc = 5000,
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

  sum_Ssq_qh_total <- matrix(0.0, nrow = 150, ncol = Q)
  n_in_cluster_qh_total <- matrix(0.0, nrow = 150, ncol = Q)
  Sum_S_qh <- matrix(0.0, nrow = 150, ncol = Q)


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

    if (round(imcmc / print_freq) == (imcmc / print_freq) ){
      print(paste0("ITER ", imcmc))
      cat("Means are:")
      print(S0_DPM$miu_h[1:actual_cluster_count])
      cat("Variances are:")
      print(S0_DPM$sigma_sq_h[1:actual_cluster_count])
      cat("Cluster counts are:")
      print(n_in_cluster[1:actual_cluster_count])
      print("IC level variance terms are:")
      print((sigma_sq_q))
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

    start_time <- Sys.time()
    # Collect updated samples from all nodes
    Aall[,] = abind::abind(clusterEvalQ(cl, A), along=1)
    elapsed_time <- Sys.time() - start_time
    timing[which(timing$Name == "mixing matrix collection"),]$Value =
      timing[which(timing$Name == "mixing matrix collection"),]$Value + as.numeric(elapsed_time)

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

    sticks_br_weights_start_time <- Sys.time()
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
      #print("NEED TO DRAW MORE STICKS!")
    }

    # Update sticks on clusters
    #export_to_cluster(cl, 'sticks', S0_DPM$sticks)
    #stick_temp <- S0_DPM$sticks
    update_cluster_vector_variable(cl, 'sticks', S0_DPM$sticks)

    elapsed_time <- Sys.time() - sticks_br_weights_start_time
    timing[which(timing$Name == "Stick Breaking Weights Time"),]$Value =
      timing[which(timing$Name == "Stick Breaking Weights Time"),]$Value + as.numeric(elapsed_time)


    # Sample U
    S0_U_start_time <- Sys.time()
    clusterEvalQ(cl, sample_u(u, sticks, cluster_memberships))
    elapsed_time <- Sys.time() - S0_U_start_time
    timing[which(timing$Name == "U Time"),]$Value =
      timing[which(timing$Name == "U Time"),]$Value + as.numeric(elapsed_time)

    # Sample cluster memberships
    cluster_membership_start_time <- Sys.time()
    export_to_cluster(cl, 'H', n_cluster)
    update_cluster_vector_variable(cl, 'n_in_cluster', n_in_cluster)
    clusterEvalQ(cl, sample_cluster_membership_indicators(cluster_memberships,
                                                          n_in_cluster,
                                                          S0, u, sticks, mu_h,
                                                          sigma_sq_h, sigma_sq_q, H))


    # Summarize cluster counts across nodes
    n_in_cluster <- Reduce(`+`, clusterEvalQ(cl, n_in_cluster))
    update_cluster_vector_variable(cl, 'n_in_cluster', n_in_cluster)
    elapsed_time <- Sys.time() - cluster_membership_start_time
    timing[which(timing$Name == "Cluster Membership Update Time"),]$Value =
      timing[which(timing$Name == "Cluster Membership Update Time"),]$Value + as.numeric(elapsed_time)

    # Cleanup
    cluster_cleanup_start_time <- Sys.time()
    generate_cluster_mapping_vector_inplace(cluster_mappings, n_in_cluster)
    export_to_cluster(cl, 'mapping', cluster_mappings)
    #update_cluster_vector_variable(cl, 'mapping', cluster_mappings)

    actual_cluster_count <- sum(cluster_mappings > 0)
    if (cluster_mappings[1] == 0){actual_cluster_count <- actual_cluster_count + 1}

    export_to_cluster(cl, 'actual_cluster_count', actual_cluster_count)

    clusterEvalQ(cl, cleanup_cluster_ordering_inplace(cluster_memberships, mu_h,
                                                      sigma_sq_h, n_in_cluster,
                                                      mapping, actual_cluster_count))

    # Cleanup - reorder the mu, sigma for the DPM
    S0_DPM$miu_h <- clusterEvalQ(cl[1], mu_h)[[1]]
    S0_DPM$sigma_sq_h <- clusterEvalQ(cl[1], sigma_sq_h)[[1]]
    n_in_cluster <- clusterEvalQ(cl[1], n_in_cluster)[[1]]
    elapsed_time <- Sys.time() - cluster_cleanup_start_time
    timing[which(timing$Name == "Cluster Cleanup Time"),]$Value =
      timing[which(timing$Name == "Cluster Cleanup Time"),]$Value + as.numeric(elapsed_time)

    elapsed_time <- Sys.time() - cluster_membership_start_time
    timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value =
      timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample S0 and Beta
    #

    spatial_map_start_time <- Sys.time()

    export_to_cluster(cl, 'mu_h_div_sigma_sq_h', S0_DPM$miu_h / S0_DPM$sigma_sq_h)

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
    suffstat_start_time <- Sys.time()
    clusterEvalQ(cl, calculate_DPM_hyperparameter_suffstats_inplace(sum_Ssq_qh,
                                                                   Sbar_qh,
                                                                   n_in_cluster_qh,
                                                                   S0,
                                                                   cluster_memberships))
    elapsed_time <- Sys.time() - suffstat_start_time
    timing[which(timing$Name == "SuffStats Calculate Time"),]$Value =
      timing[which(timing$Name == "SuffStats Calculate Time"),]$Value + as.numeric(elapsed_time)

    # Unite the sufficient statistics across the different nodes
    suffstat_pull_cluster_start_time <- Sys.time()
    utility_update_matrix_across_nodes_inplace(sum_Ssq_qh_total, clusterEvalQ(cl, sum_Ssq_qh))
    utility_update_matrix_across_nodes_inplace(n_in_cluster_qh_total, clusterEvalQ(cl, n_in_cluster_qh))
    utility_update_matrix_across_nodes_inplace(Sum_S_qh, clusterEvalQ(cl, Sbar_qh))

    # Account for chance that none in cluster -> can cause nan
    n_in_cluster_qh_total_temp <- n_in_cluster_qh_total
    n_in_cluster_qh_total_temp[n_in_cluster_qh_total_temp == 0] <- 1
    #Sum_S_qh              <- Reduce(`+`, clusterEvalQ(cl, Sbar_qh) )
    Sbar_qh_total         <- Sum_S_qh / n_in_cluster_qh_total_temp # also not scaled by variance yet

    # Now obtain the S - Sbar squared term - specific to each node
    # TODO - get rid of this quantity
    DPM_hyper_sse_unscaled <- (Sum_S_qh^2 - 2*Sum_S_qh*Sbar_qh_total + Sbar_qh_total^2)

    elapsed_time <- Sys.time() - suffstat_pull_cluster_start_time
    timing[which(timing$Name == "SuffStats Node Summary Time"),]$Value =
      timing[which(timing$Name == "SuffStats Node Summary Time"),]$Value + as.numeric(elapsed_time)

    hyperparameter_sample_start_time <- Sys.time()
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

    # Push new cluster parameters to each worker process
    update_cluster_vector_variable(cl, 'mu_h', S0_DPM$miu_h)
    update_cluster_vector_variable(cl, 'sigma_sq_h', S0_DPM$sigma_sq_h)
    elapsed_time <- Sys.time() - hyperparameter_sample_start_time
    timing[which(timing$Name == "Sample cluster parameters Time"),]$Value =
      timing[which(timing$Name == "Sample cluster parameters Time"),]$Value + as.numeric(elapsed_time)



    if (any(is.na(S0_DPM$miu_h))){
      stop("NA cluster mean obtained")
    }

    # Update concentration parameter
    concentration_start_time <- Sys.time()
    phi <- rbeta(1, concentration + 1, Q*V)
    odds <- (4 + H - 1) / (Q*V * (4.0 - log(phi)))

    if (runif(1, 0,1) < (odds/(odds+1)) ){
      concentration <- rgamma(1, 4.0 + H, scale = 1.0 / (4.0 - log(phi)))
    } else {
      concentration <- rgamma(1, 4.0 + H - 1.0, scale = 1.0 / (4.0 - log(phi)))
    }

    elapsed_time <- Sys.time() - concentration_start_time
    timing[which(timing$Name == "Sample concentration time"),]$Value =
      timing[which(timing$Name == "Sample concentration time"),]$Value + as.numeric(elapsed_time)


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
    sse_start_time <- Sys.time()
    clusterEvalQ(cl, evaluate_sse_inplace(sse, Ei))
    elapsed_time <- Sys.time() - sse_start_time
    timing[which(timing$Name == "SSE Time"),]$Value =
      timing[which(timing$Name == "SSE Time"),]$Value + as.numeric(elapsed_time)


    # Evaluate the regression scale term
    reg_scale_time <- Sys.time()
    clusterEvalQ(cl, evaluate_reg_scale_inplace(reg_scale, S0, Beta, mu_h, sigma_sq_h,
                                               tau_sq, cluster_memberships, lambda_sq))
    elapsed_time <- Sys.time() - reg_scale_time
    timing[which(timing$Name == "Scale Time"),]$Value =
      timing[which(timing$Name == "Scale Time"),]$Value + as.numeric(elapsed_time)

    # Aggregate across worker processes
    cluster_process_reduction_start_time <- Sys.time()
    sse_total       <- Reduce(`+`, clusterEvalQ(cl, sse) )
    reg_scale_total <- Reduce(`+`, clusterEvalQ(cl, reg_scale) )
    elapsed_time <- Sys.time() - cluster_process_reduction_start_time
    timing[which(timing$Name == "Sigma2 Node Reduction Time"),]$Value =
      timing[which(timing$Name == "Sigma2 Node Reduction Time"),]$Value + as.numeric(elapsed_time)

    # Draw new variance parameters
    sigma_sq_update_time <- Sys.time()
    shape      = (N * V*Q + Q*(P+1) * V) / 2.0
    scale      = sum(sse_total) / 2.0 + sum(reg_scale_total) / 2.0
    sigma_sq_q[1:Q] = 1 / rgamma(Q, shape, rate = scale)
    export_to_cluster(cl, 'sigma_sq_q', sigma_sq_q)
    elapsed_time <- Sys.time() - sigma_sq_update_time
    timing[which(timing$Name == "Sigma Sq Q"),]$Value =
      timing[which(timing$Name == "Sigma Sq Q"),]$Value + as.numeric(elapsed_time)

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
    local_shrinkage_start_time <- Sys.time()
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
    elapsed_time <- Sys.time() - local_shrinkage_start_time
    timing[which(timing$Name == "Local Shrinkage"),]$Value =
      timing[which(timing$Name == "Local Shrinkage"),]$Value + as.numeric(elapsed_time)


    # Global (Group) Shrinkage Parameters
    global_shrinkage_start_time <- Sys.time()
    clusterEvalQ(cl, calculate_tausum_inplace(tau_sum, Beta, lambda_sq, sigma_sq_q))

    tau_sum[,] = Reduce(`+`, clusterEvalQ(cl, tau_sum) )

    # Specific Tau (group horseshoe)
    # for (q in 1:Q){
    #   for (p in 1:P){
    #     tau_sq[q, p] = max(invgamma::rinvgamma(1,  (V+1)/2, rate = 1/xi[q,p] + tau_sum[q,p] ), 1e-20)
    #     xi[q, p]     =  invgamma::rinvgamma(1,  1, rate = 1 + 1/tau_sq[q,p] )
    #   }
    # }
    # Common Tau
    tau_sq[,] = max(invgamma::rinvgamma(1,  (V*P*Q+1)/2, rate = 1/sum(xi[1,1]) + sum(tau_sum) ), 1e-20)
    xi[,]     = invgamma::rinvgamma(1,  1, rate = 1 + 1/tau_sq[1,] )


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
    elapsed_time <- Sys.time() - local_shrinkage_start_time
    timing[which(timing$Name == "Global Shrinkage"),]$Value =
      timing[which(timing$Name == "Global Shrinkage"),]$Value + as.numeric(elapsed_time)



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
    Atemp <- Ai_posterior_mean[ ((i-1)*Q+1):(i*Q), ]
    Ai_3D[,,i] <- Atemp %*% Atemp %*% sqrtm(solve(t(Atemp)%*%Atemp));
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
    Beta_pseudo_pvalue = Beta_pvalues_reordered,
    sigma_sq_q_posterior_mean = sigma_sq_q_posterior_mean,
    tau_sq_posterior_mean = tau_sq_posterior_mean,
    lambda_sq_posterior_mean = lambda_sq_PM_reordered,
    timing = timing)
  )

}

