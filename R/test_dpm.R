if (0 == 1){


  #
  # Version 1: Hyperparameter Check ----
  #


  library(hcicaR)

  # data path
  dp <- "/Users/joshlukemire/Documents/Research/cbis/SparseBayesICA/simulation_demo_data"

  set.seed(100)

  N <- 10
  V <- 100000
  Q <- 3
  Htrue <- 5
  mu_h_true       <- c(2.0, 0.0, -3.0, 0.5, -0.5)
  sigma_h_sq_true <- c(2.0, 0.1, 0.5, 0.25, 0.8)
  CMtrue <- matrix(as.integer(sample(1:Htrue, size = Q*V, replace = TRUE, prob = c(100, 100, 30, 10, 5) )),
                   nrow = V,
                   ncol  = Q)

  # Generate true values for S0
  S0true <- matrix(0, nrow = V, ncol = Q)
  for (q in 1:Q){
    for (v in 1:V){
      S0true[v, q] <- rnorm(1, mean = mu_h_true[CMtrue[v, q]], sd = sqrt(sigma_h_sq_true[CMtrue[v, q]]) )
    }
  }

  # Setup intermediate quantities
  sum_Ssq_qh      <- matrix(0, nrow = 50, ncol = Q)
  Sbar_qh         <- matrix(0, nrow = 50, ncol = Q)
  n_in_cluster_qh     <- matrix((0), nrow = 50, ncol = Q)
  S0                  <- S0true + 0
  cluster_memberships <- matrix(as.integer(CMtrue - as.integer(1)), nrow = V)

  mu_h <- rep(0.0, Htrue)
  sigma_sq_h <- rep(1.0, Htrue)

  n_in_cluster <- rep(0, Htrue)
  for (h in 1:Htrue){
    n_in_cluster[h] <- sum(CMtrue == h)
  }

  ##########################################################
  #

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

  prior_shape = 1.0
  prior_scale = 1.0
  sigma_sq_q <- rep(0.7, Q)

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
                                     Htrue)

  ###
  ### Normal inverse gamma version
  # ###
  # s0_qh_bar <- matrix(0, nrow = Htrue, ncol = Q)
  # for (h in 1:Htrue){
  #   for (q in 1:Q){
  #     s0_qh_bar[h, q] <- var(S0[cluster_memberships[,q] == (h-1), q ]) * (sum(sigma_sq_q) /  sigma_sq_q[q])
  #   }
  # }

  print(data.frame(mu = mu_h, sigma_sq = sigma_sq_h ))


  #
  # Version 2: Cluster Assignments Check ----
  #

  set.seed(100)

  N <- 10
  V <- 1000
  Q <- 3
  Htrue <- 5
  mu_h_true       <- c(12.0, 0.0, -13.0, 5.5, -5.5)
  sigma_h_sq_true <- c(2.0, 0.1, 0.5, 0.25, 0.8)
  CMtrue <- matrix(as.integer(sample(1:Htrue, size = Q*V, replace = TRUE, prob = c(100, 100, 30, 10, 5) )),
                   nrow = V,
                   ncol  = Q)

  # Generate true values for S0
  S0true <- matrix(0, nrow = V, ncol = Q)
  for (q in 1:Q){
    for (v in 1:V){
      S0true[v, q] <- rnorm(1, mean = mu_h_true[CMtrue[v, q]], sd = sqrt(sigma_h_sq_true[CMtrue[v, q]]) )
    }
  }

  #### Initial Values for Test
  cluster_memberships <- matrix(as.integer(sample(0:10, size = Q*V, replace = TRUE)), nrow = V)
  sticks <- cumprod(runif(50))
  u <- matrix(runif(Q*V, min = 0, max = sticks[cluster_memberships]), nrow = V)
  concentration = 1.0

  n_in_cluster <- rep(0, 50)
  for (h in 1:50){
    n_in_cluster[h] <- sum(cluster_memberships == h)
  }

  mu_h       <- c(mu_h_true, rnorm(45))
  sigma_sq_h <- c(sigma_h_sq_true, runif(45))

  S0 <- S0true + 0.0

  cluster_mappings = rep(0, length(mu_h))


  ###################

  for (j in 1:10000){

    print(j)

  # Get the minimum value of u (across all nodes)
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

  print(max(cluster_memberships))

  }

  table(CMtrue, cluster_memberships)




  #
  # Version 3: S0 Estimation Accuracy Check ----
  #

  sample_spatial_maps_inplace(S0,
                              Beta,
                              posterior_mean,
                              posterior_variance,
                              prior_precision,
                              draw,
                              AtY,
                              X_with_int,
                              XtX,
                              sigma_sq_q,
                              tau_sq,
                              lambda_sq,
                              mu_h_div_sigma_sq_h,
                              sigma_sq_h,
                              cluster_memberships,
                              Ip)




}
