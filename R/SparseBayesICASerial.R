#---- Main Function
# Y <- Q x V x N
# initial_values should contain A, S0, and beta
SparseBayesICASerial <- function(Y, X, initial_values,
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
  angle_range = seq(from = 0.0, to = 2*pi, length.out = 500)
  sin_theta <- sin(angle_range)
  cos_theta <- cos(angle_range)

  # Storage for suff stat for mixing matrix update
  QxQStore    <- matrix(0, nrow = Q, ncol = Q)
  VxQStore    <- matrix(0, nrow = V, ncol = Q)
  YSiall      <- array(0, dim = c(Q*N, Q))
  YSiOmega    <- array(0, dim = c(Q*N, Q))
  Aall        <- initial_values$A + 0
  AtY         <- initial_values$Si + 0
  X_with_int  <- cbind(1, X)
  XtX         <- t(X_with_int) %*% X_with_int
  Ip <- diag(1, ncol(X_with_int))
  posterior_mean     <- rep(0, ncol(X_with_int))
  posterior_variance <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  prior_precision    <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  draw <- rep(0, ncol(X_with_int))
  sse <- rep(0, Q)
  reg_scale <- rep(0, Q)
  column_reduced_mixing_matrix <- matrix(0, nrow = Q, ncol = Q - 2)
  null_space <- matrix(0, nrow = Q, ncol = 2)

  YYt <- matrix(0, nrow=N*Q, Q)
  for (iSubj in 1:N){
    inds <-  (Q*(iSubj-1)+1):(Q*iSubj)
    YYt[ inds, ] = t(Y[,inds]) %*% Y[,inds]
  }

  S0_DPM = initial_values$S0_DPM
  cluster_mappings = rep(0, length(S0_DPM$n_in_cluster))
  n_in_cluster <- S0_DPM$n_in_cluster + 0

  cluster_memberships <- S0_DPM$cluster_membership
  if (min(cluster_memberships) == 1){
    cluster_memberships = cluster_memberships - as.integer(1)
  }

  u          = S0_DPM$u
  sticks     = S0_DPM$sticks
  mu_h       = S0_DPM$miu_h
  sigma_sq_h = S0_DPM$sigma_sq_h

  sum_Ssq_qh      <- matrix(0, nrow = 50, ncol = Q)
  Sbar_qh         <- matrix(0, nrow = 50, ncol = Q)
  n_in_cluster_qh <- matrix(0, nrow = 50, ncol = Q)

  # Parameters
  S0         <- initial_values$S0 + 0
  Si         <- initial_values$Si + 0
  sigma_sq_q <- initial_values$sigma_sq_q + 0
  Beta       <- initial_values$beta + 0
  tau_sum    <- matrix(0, nrow = Q, ncol = P)
  tau_sq     <- initial_values$Horseshoe$tau_sq + 0
  xi         <- initial_values$Horseshoe$xi + 0
  lambda_sq  <- initial_values$Horseshoe$lambda_sq + 0
  Ei         <- matrix(0, nrow = V, ncol = Q * N)
  nu         <- matrix(0.1, nrow = V, ncol = Q*P)
  cauchy_mixing       <- matrix(0.1, nrow = V, ncol = Q*P)
  cauchy_mixing_prior <- matrix(0.1, nrow = V, ncol = Q*P)

  sigma_sq_q = initial_values$sigma_sq_q + 0

  # Storage - posterior quantities
  Ai_posterior_mean <- array(0, dim = c(Q*N, Q))
  sigma_sq_q_posterior_mean <- rep(0, Q)
  tau_sq_posterior_mean <- matrix(0, nrow = Q, ncol = P)
  S0_posterior_mean  <- matrix(0, nrow = V, ncol = Q)
  Beta_posterior_mean <- matrix(0, nrow = V, ncol = Q*P)
  Beta_ngt0 <- matrix(0, nrow = V, ncol = Q*P)
  lambda_sq_posterior_mean <- matrix(0, nrow = V, ncol = Q*P)

  # Verify input dimensions all match

  # Initialize object to keep track of computation time
  timing <- create_sparsebayes_timing_log()

  # Obtain the cluster splits

  overall_start_time <- Sys.time()
  for (imcmc in 1:(nburnin + nmcmc)){

    if (round(imcmc / print_freq) == (imcmc / print_freq) ){
      print(paste0("ITER ", imcmc))
    }

    # Store prior precision
    Omega <- diag(1.0 / sigma_sq_q)

    #
    # Mixing Matrix Update
    #

    # Monitor Mixing Matrix Update Time
    mixing_matrix_start_time = Sys.time()

    #ttt <- t(Y[,1:5]) %*% Si[,1:5]

    #Get the required quantities for the mbvmf update
    calculate_YiSit_inplace(YSiall, QxQStore, Y, Si, Q)

    # print("Si is it updated???")
    # print(YSiall[1:5, 1:5])

    YSiOmega[,]    <- YSiall %*% Omega

    if (any(is.na(YSiOmega)) || any(is.na(YSiall))){
      stop("Missing values in Y x Si")
    }

    endIndex = 0
    for (i in 1:N){
      startIndex = endIndex + 1
      endIndex   = endIndex + Q
      negYSiSigInvover2 <- -Omega/2
      Aall[startIndex:endIndex,] <- gibbs_sample_mixing_matrix(Aall[startIndex:endIndex,],
                                                            (YSiall[startIndex:endIndex,]) %*% Omega,
                                                            (negYSiSigInvover2),
                                                            YYt[startIndex:endIndex,],
                                                            angle_range, sin_theta,
                                                            cos_theta,
                                                            column_reduced_mixing_matrix,
                                                            null_space)
    }


    # Finish Monitoring Mixing Matrix Update Time
    elapsed_time <- Sys.time() - mixing_matrix_start_time
    timing[which(timing$Name == "Overall Mixing Matrix Time"),]$Value =
      timing[which(timing$Name == "Overall Mixing Matrix Time"),]$Value + as.numeric(elapsed_time)

    #
    # AtY Bookkeeping - is Ai' Yi
    #

    AtY_start_time <- Sys.time()

    # print("in")
    # print(AtY[1:3, 1:3])
    calculate_AitYi_inplace(AtY, VxQStore, Aall, Y, Q)
    # print("out")
    # print(AtY[1:3, 1:3])

    elapsed_time <- Sys.time() - AtY_start_time
    timing[which(timing$Name == "At x Y Time"),]$Value =
      timing[which(timing$Name == "At x Y Time"),]$Value + as.numeric(elapsed_time)
    timing[which(timing$Name == "Overall Variable Management Time"),]$Value =
      timing[which(timing$Name == "Overall Variable Management Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample DPM Quantities
    #

    cluster_membership_start_time <- Sys.time()

    # Get the minimum value of u
    min_u =  min(u)

    # Sample stick breaking weights
    H <- which(n_in_cluster == 0)[1] - 1
    n_cluster      <- H + 10
    n_remaining    <- Q*V - cumsum(n_in_cluster[1:n_cluster])
    draws          <- rbeta(n_cluster, 1 + n_in_cluster[1:n_cluster], concentration + n_remaining)
    stick_remaining <- cumprod(c(1, 1 - draws))
    mu_h[(H+1):n_cluster]      <- rnorm(10)
    sigma_sq_h[(H+1):n_cluster] <- rinvgamma(10, 1, 1)
    sticks[1:n_cluster] <- draws * stick_remaining[-length(draws)]

    # Draw any additional required clusters
    if (min_u < sticks[n_cluster]){
      #print("Min sticks was:")
      #print(S0_DPM$sticks[n_cluster])
      #print("NEED TO DRAW MORE STICKS!")
    }

    # Sample U
    sample_u(u, sticks, cluster_memberships)

    # Sample cluster memberships
    sample_cluster_membership_indicators(cluster_memberships,
                                         n_in_cluster,
                                         S0, u, sticks, mu_h,
                                         sigma_sq_h, sigma_sq_q, n_cluster)


    # Cleanup
    generate_cluster_mapping_vector_inplace(cluster_mappings, n_in_cluster)

    actual_cluster_count <- sum(cluster_mappings > 0)
    if (cluster_mappings[1] == 0){actual_cluster_count <- actual_cluster_count + 1}

    cleanup_cluster_ordering_inplace(cluster_memberships, mu_h,
                                     sigma_sq_h, n_in_cluster,
                                     cluster_mappings, actual_cluster_count)

    elapsed_time <- Sys.time() - cluster_membership_start_time
    timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value =
      timing[which(timing$Name == "Overall Cluster Membership Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample S0 and Beta
    #

    spatial_map_start_time <- Sys.time()

    mu_h_div_sigma_sq_h = mu_h / sigma_sq_h

    sample_spatial_maps_inplace(S0, Beta, posterior_mean,
                                posterior_variance,
                                prior_precision, draw,
                                AtY, X_with_int, XtX, sigma_sq_q,
                                tau_sq, lambda_sq,
                                mu_h_div_sigma_sq_h, sigma_sq_h,
                                cluster_memberships,
                                Ip)

    if (0 == 1){
      ts <- seq(from = 1, to = 125, by = 5)
      temp <- AtY[,ts]
      s0hat <- apply(temp, 1, mean)
      plot(s0hat, S0[,1])
    }

    if (0 == 1){
      for (q in 1:Q){
        ts <- seq(from = q, to = 125, by = Q)
        temp <- AtY[,ts]
        s0hat <- apply(temp, 1, mean)
        S0[,q] <- s0hat
        Beta = Beta * 0.00000000000000000001
      }
      print(S0[1:3, 1:4])
    }


    if ( any(is.na(Beta)) ){
     stop("NA beta was obtained")
    }


  elapsed_time <- Sys.time() - spatial_map_start_time
  timing[which(timing$Name == "Overall Spatial Map Time"),]$Value =
    timing[which(timing$Name == "Overall Spatial Map Time"),]$Value + as.numeric(elapsed_time)

    #
    # Sample DPM Hyperparameters
    #

    DPM_hyperparameter_time <- Sys.time()

    # Evaluate the "sufficient statistics"
    calculate_DPM_hyperparameter_suffstats_inplace(sum_Ssq_qh,
                                                   Sbar_qh,
                                                   n_in_cluster_qh,
                                                   S0,
                                                   cluster_memberships)

    # Unite the sufficient statistics across the different nodes
    sum_Ssq_qh_total      <- sum_Ssq_qh
    n_in_cluster_qh_total <- n_in_cluster_qh

    # Account for chance that none in cluster -> can cause nan
    n_in_cluster_qh_total_temp <- n_in_cluster_qh_total
    n_in_cluster_qh_total_temp[n_in_cluster_qh_total_temp == 0] <- 1
    Sum_S_qh              <- Sbar_qh
    Sbar_qh_total         <- Sum_S_qh / n_in_cluster_qh_total_temp # also not scaled by variance yet

    # Now obtain the S - Sbar squared term - specific to each node
    # TODO - get rid of this quantity
    DPM_hyper_sse_unscaled <- (Sum_S_qh^2 - 2*Sum_S_qh*Sbar_qh_total + Sbar_qh_total^2)

    sample_DPM_hyperparameters_inplace(mu_h,
                                       sigma_sq_h,
                                       sum_Ssq_qh_total,
                                       Sbar_qh_total,
                                       DPM_hyper_sse_unscaled,
                                       n_in_cluster_qh_total,
                                       n_in_cluster,
                                       sigma_sq_q,
                                       prior_shape,
                                       prior_scale,
                                       actual_cluster_count)

    # cat("Means are:")
    # print(mu_h[1:actual_cluster_count])
    # cat("Variances are are:")
    # print(sigma_sq_h[1:actual_cluster_count])
    # cat("Cluster counts are:")
    # print(n_in_cluster[1:actual_cluster_count])

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
    evaluate_error_inplace(Ei, Si, AtY, S0, Beta, X_with_int)

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
    evaluate_sse_inplace(sse, Ei)
    # cat("SSE Ei")
    # print(sse)
    # evaluate_ssY_inplace(sse, AtY)
    # print("AtY after SSY function")
    # print(AtY[1:3, 1:3])
    # cat("SSE AY")
    # print(sse)

    # Evaluate the regression scale term
    evaluate_reg_scale_inplace(reg_scale, S0, Beta, mu_h, sigma_sq_h,
                               tau_sq, cluster_memberships, lambda_sq)




    # Aggregate across worker processes
    sse_total       <- sse
    reg_scale_total <- reg_scale

    # Draw new variance parameters
    # shape      = (N * V + (P+1) * V) / 2.0
    # scale      = sse_total / 2.0 + reg_scale_total / 2.0
    # sigma_sq_q = 1 / rgamma(Q, shape, rate = scale)
    shape      = (N * V*Q + Q*(P+1) * V) / 2.0
    scale      = sum(sse_total) / 2.0 + sum(reg_scale_total) / 2.0
    sigma_sq_q[1:Q] = 1 / rgamma(Q, shape, rate = scale)



    #print("Draw is:")
    #print(sigma_sq_q)



    # Check if failed
    if (any(is.na(sigma_sq_q))){
      print("Sigma squared h")
      print(sigma_sq_h)
      print("n in cluster")
      print(n_in_cluster)
      print("SSE")
      print(sse_total)
      print("Reg Scale term:")
      print(reg_scale_total)
      stop("NA variance obtained")
    }

    # Local Shrinkage Parameters
    sample_local_shrinkage_parameters_inplace(lambda_sq,
                                              nu,
                                              cauchy_mixing,
                                              cauchy_mixing_prior,
                                              Beta,
                                              tau_sq,
                                              sigma_sq_q)

    if (any(is.na(lambda_sq))){
      stop("NA Lambda obtained")
    }

    # Global (Group) Shrinkage Parameters
    calculate_tausum_inplace(tau_sum, Beta, lambda_sq, sigma_sq_q)

    for (q in 1:Q){
      for (p in 1:P){
        tau_sq[q, p] = max(1 / rgamma(1,  (V+1)/2, rate = 1/xi[q,p] + tau_sum[q,p] ), 1e-200)
        xi[q, p]     = 1 / rgamma(1,  1, rate = 1 + 1/tau_sq[q,p] )
      }
    }

    if (any(is.na(tau_sq))){
      stop("NA Tau-squared obtained")
    }

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

      update_spatial_map_posterior_tracking_inplace(S0_posterior_mean,
                                                    Beta_posterior_mean,
                                                    Beta_ngt0,
                                                    lambda_sq_posterior_mean,
                                                    S0, Beta, lambda_sq, nmcmc)



      elapsed_time <- Sys.time() - sample_storage_time
      timing[which(timing$Name == "Overall Sample Storage Time"),]$Value =
        timing[which(timing$Name == "Overall Sample Storage Time"),]$Value + as.numeric(elapsed_time)



    }

  }

  overall_run_time <- Sys.time() - overall_start_time
  timing[which(timing$Name == "Overall MCMC Time"),]$Value = as.numeric(overall_run_time)

  # Reshape variables to more easily understood format
  Ai_3D <- array(0, dim = c(Q, Q, N))
  for (i in 1:N){
    Ai_3D[,,i] <- Ai_posterior_mean[ ((i-1)*Q+1):(i*Q), ]
  }

  S0_PM   <- S0_posterior_mean

  # Get p-values based on credible intervals
  Beta_prop_gt0 <- Beta_ngt0 / nmcmc
  Beta_pvals <- 1 - Beta_prop_gt0
  Beta_pvals[Beta_pvals > 0.5] <- 1.0 - Beta_pvals[Beta_pvals > 0.5]
  Beta_pvals <- 2.0 * Beta_pvals # two-sided correction

  # Reorder Betas and P-values into 3d array, where 3rd dimension indexes covariates
  Beta_PM_reordered <- array(0, dim = c(V, Q, P))
  lambda_sq_PM_reordered <- array(0, dim = c(V, Q, P))
  Beta_pvalues_reordered <- array(0, dim = c(V, Q, P))
  for (q in 1:Q){
    index_set <- ((q-1)*P+1):(q*P)
    Beta_PM_reordered[,q,]      <- Beta_posterior_mean[, index_set]
    Beta_pvalues_reordered[,q,] <- Beta_pvals[, index_set]
    lambda_sq_PM_reordered[,q,] <- lambda_sq_posterior_mean[, index_set]
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

