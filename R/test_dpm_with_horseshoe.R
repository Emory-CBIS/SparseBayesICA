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

  P <- 2
  tau_sq_true    <- matrix(c(0.01, 0.001), nrow = Q, ncol = P)
  lambda_sq_true <- matrix(sample(c(0.01, 1000), size = Q * V * P, replace = TRUE, prob = c(0.9, 0.1)), nrow = V)
  # Generate the betas
  beta_true <- matrix(0, nrow = V, ncol = Q*P)
  cind <- 0
  for (q in 1:Q){
    for (p in 1:P){
      cind = cind + 1
      for (v in 1:V){
        beta_true[v,cind] <- rnorm(1, mean = 0, sd = sqrt(tau_sq_true[q,p] * lambda_sq_true[v, cind]) )
      }
    }
  }

  X <- matrix(rnorm(N*P), nrow = N)
  X_with_int  <- cbind(1, X)
  XtX         <- t(X_with_int) %*% X_with_int
  posterior_mean     <- rep(0, ncol(X_with_int))
  posterior_variance <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  prior_precision    <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  draw <- rep(0, ncol(X_with_int))

  # Generate the "observed" data, this is directly A'Y
  AtY <- matrix(0, nrow = V, ncol = N * Q)
  inds <- 1:Q
  for (i in 1:N){
    AtYi <- S0true
    for (p in 1:P){
      bs <- seq(from = p, to = ncol(beta_true), by = P)
      AtYi <- AtYi + X[i, p] * beta_true[, bs]
    }
    AtY[,inds] <- AtYi
    inds <- inds + Q
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

  # input quantities
  S0   <- matrix(0, nrow = V, ncol = Q)
  Beta <- matrix(0, nrow = V, ncol = P*Q)
  tau_sq <- tau_sq_true
  lambda_sq <- lambda_sq_true + 0
  cluster_memberships <- CMtrue - as.integer(1)
  mu_h_div_sigma_sq_h <- mu_h_true / sigma_h_sq_true
  sigma_sq_h = sigma_h_sq_true
  Ip <- diag(1, ncol(X_with_int))


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

  plot(S0[,1], S0true[,1])
  plot(S0[,2], S0true[,2])
  plot(S0[,3], S0true[,3])

  plot(Beta[,1], beta_true[,1])
  plot(Beta[,5], beta_true[,5])





  ###### now hyperparameters

  nu = 1 / lambda_sq_true
  cauchy_mixing = lambda_sq_true + 0
  cauchy_mixing_prior = 1 / lambda_sq_true

  # Local Shrinkage Parameters
  sample_local_shrinkage_parameters_inplace(lambda_sq,
                                            nu,
                                            cauchy_mixing,
                                            cauchy_mixing_prior,
                                            beta_true,
                                            tau_sq_true,
                                            sigma_sq_q/sigma_sq_q)


  plot(lambda_sq[,1], lambda_sq_true[,1])

  ttt <- lambda_sq[lambda_sq_true == 1000]
  median(ttt)

  if (any(is.na(lambda_sq))){
    stop("NA Lambda obtained")
  }

  # Global (Group) Shrinkage Parameters
  tau_sum  <- matrix(0, nrow = Q, ncol = P)
  tau_sq <- matrix(0, nrow = Q, ncol = P)
  xi <- matrix(0.1, nrow = Q, ncol = P)

  calculate_tausum_inplace(tau_sum, beta_true, lambda_sq_true, sigma_sq_q)

  for (q in 1:Q){
    for (p in 1:P){
      tau_sq[q, p] = max(1 / rgamma(1,  (V+1)/2, rate = 1/xi[q,p] + tau_sum[q,p] ), 1e-200)
      xi[q, p]     = 1 / rgamma(1,  1, rate = 1 + 1/tau_sq[q,p] )
    }
  }

  tau_sq
  tau_sq_true









}
