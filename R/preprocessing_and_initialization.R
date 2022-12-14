# Codes Related to Initial Decomposition and Guess

# Function to obtain prewhitened data
# TODO, make sure pipe works + output is formatted correctly (ideally take nifti path as an argument)
preprocess_subject <- function(data, Q, unit_variance = FALSE, center = FALSE, include_sigma_ml = TRUE){

  # Dimensions
  Ti = dim(data)[1]
  V  = dim(data)[2]

  # Center voxels, this is done in case previous processing gave them some other mean
  # (for example, mean = 100 is common output from some pipelines)
  data <- sweep(data, 2, colMeans(data))

  # Center time points
  data <- sweep(data, 1, rowMeans(data))

  # PCA dimension reduction
  data_pca <- prcomp(t(data), center = center)

  # Get the eigenvals
  eigvals_Q <- data_pca$sdev[1:Q]^2
  mean_remaining_eigvals <- mean(data_pca$sdev[-c(1:Q)]^2)
  #sum_remaining_eigvals <- sum(data_pca$sdev[-c(1:Q)]^2)

  # Eigenvectors
  U_eigen <- data_pca$rotation[, 1:Q]

  # Construct the transformation matrix
  tsfm_matrix <- t(U_eigen)
  if (unit_variance == TRUE){
    #tsfm_matrix <- diag( (eigvals_Q - mean_remaining_eigvals)^(-1/2) ) %*% tsfm_matrix
    tsfm_matrix <- diag( (eigvals_Q)^(-1/2) ) %*% tsfm_matrix
    if (include_sigma_ml == TRUE){
      tsfm_matrix <- diag(rep(1 + mean_remaining_eigvals, Q))^(1/2) %*% tsfm_matrix
      #tsfm_matrix <- diag(rep(1 + sum_remaining_eigvals, Q))^(1/2) %*% tsfm_matrix
    }
  }


  # Construct the whitened data and flip it to be V x Q
  Y = t(tsfm_matrix %*% data)

  # Make sure the whitened data are centered (time courses mean 0)
  Y <- sweep(Y, 2, colMeans(Y))


  return_list <- list(Y = Y,
                      eigenvalues = data_pca$sdev^2)

  # Returns whitened data
  return(return_list)

}


PCA_dimension_reduction <- function(data, nPC){

  # Dimensions
  Ti = dim(data)[1]
  V  = dim(data)[2]

  # Center time points
  data <- sweep(data, 1, rowMeans(data))

  # PCA dimension reduction
  data_svd <- svd(t(data), nu = nPC)

  data_pc <- t(data_svd$u %*% diag(data_svd$d[1:nPC]))

  # Returns first nPC principal components
  return(data_pc)

}

obtain_s0_guess_fastica <- function(data){

  Q = dim(data)[1]
  V = dim(data)[2]

  # Apply fastica
  ica_decomposition <- fastICA(t(data), Q)

  S0 = ica_decomposition$S

  return(S0)

}



obtain_initial_guess_crosssectional <- function (data, S0, X,
                                                 nMoG = 2){

  # Dimensions
  N = nrow(X)
  V = dim(S0)[1]
  Q = dim(S0)[2]
  P = ncol(X)

  # First project prewhitened time courses onto S0 to obtain initial mixing matrix
  proj <- t(S0) %*% solve(  S0 %*% t(S0) )

  Si = array(0, dim=c(V * N*Q) )
  A = array(0, dim=c(N, Q, Q) )

  sse_level1 <- 0.0

  for (i in 1:N){
    # Project
    Atemp = data[[i]] %*% proj
    # Apply orthonormal tsfm
    A[i,,] = Atemp %*% sqrtm(solve(t(Atemp)%*%Atemp));
    # Back reconstruct Si using new A
    Si[i,,] = t(A[i,,]) %*% data[[i]]

    # Error at first level of model
    errs_i = data[[i]] - A[i,,] %*% Si[i,,]
    sse_level1 <- sse_level1 + sum(errs_i^2)
  }

  sigma_1_sq <- sse_level1 / (N * V * Q - 1)

  # Project subject-specific maps onto design matrix to obtain initial values for
  # S0 and the covariate effects
  Xaug <- cbind(1, X)
  proj = solve(t(Xaug) %*% Xaug) %*% t(Xaug)

  S0   <- matrix(0, nrow = Q, ncol = V)
  Beta <- array(0, dim = c(P, Q, V))

  sigma_2_sq <- rep(0, Q)

  for (iIC in 1:Q){
    est = proj %*% Si[,iIC,]
    S0[iIC,] = est[iIC,]
    Beta[,iIC,] = est[2:P+1,]

    errs = Si[,iIC,] - Xaug %*% est
    sigma_2_sq[iIC] = var(c(errs))

  }

  # MoG Terms
  miu_z = matrix(0, nrow = Q, ncol = nMoG)
  sigma_z_sq = matrix(0, nrow = Q, ncol = nMoG)
  pi_z = matrix(0, nrow = Q, ncol = nMoG)


  for (iIC in 1:Q){
    S0_q = S0[iIC,]
    sigma_noise = sqrt(var(S0_q))

    # Determine sign of activation mean
    quants = quantile(S0_q, c(0.025, 0.975));
    MoG_sign = 1;
    if (abs(quants[1]) > abs(quants[2])){
      MoG_sign = -1
    }

    cutpoint = 1.64 * sigma_noise
    if (sum((S0_q*MoG_sign) > cutpoint) > 0){
      if (abs(quants[1]) > abs(quants[2])){
        cutpoint = quants[1];
      } else {
        cutpoint = quants[2];
      }
    }

    miu_z[iIC, 2]      = mean( S0_q[ (S0_q*MoG_sign) > (MoG_sign*cutpoint) ] );
    sigma_z_sq[iIC, 2] = var( S0_q[ (S0_q*MoG_sign) > (MoG_sign*cutpoint) ] );
    sigma_z_sq[iIC, 1] = sigma_noise^2;
    pi_z[iIC, 2]       = sum( (S0_q*MoG_sign) > (MoG_sign*cutpoint) ) / length(S0_q);
    pi_z[iIC, 1]       = 1 - pi_z[iIC, 2];

  }

  guess = list(
    S0 = S0,
    Beta = Beta,
    A = A,
    sigma_1_sq = sigma_1_sq,
    sigma_2_sq = sigma_2_sq,
    miu_z = miu_z,
    sigma_z_sq = sigma_z_sq,
    pi_z = pi_z
  )

  return(guess)

}












obtain_initial_guess_sparsebayes <- function (data, S0, X, n_initial_cluster = 10){

  # Dimensions
  N = nrow(X)
  V = dim(S0)[1]
  Q = dim(S0)[2]
  P = ncol(X)

  # First project prewhitened time courses onto S0 to obtain initial mixing matrix
  proj <- solve(  t(S0) %*% S0 ) %*% t(S0)

  Si = array(0, dim=c(V, N*Q) )
  A  = array(0, dim=c(N*Q, Q) )

  sse_level1 <- 0.0

  eI = 0
  for (i in 1:N){

    # Project
    Atemp = t(proj %*% data[[i]])

    # Apply orthonormal tsfm
    sI = eI + 1
    eI = eI + Q
    A[sI:eI,] = Atemp %*% sqrtm(solve(t(Atemp)%*%Atemp));

    # Back reconstruct Si using new A
    Si[,sI:eI] = data[[i]] %*% A[sI:eI,]

    # Error at first level of model
    errs_i = data[[i]] - Si[,sI:eI] %*% t(A[sI:eI,])
    sse_level1 <- sse_level1 + sum(errs_i^2)
  }

  sigma_1_sq <- sse_level1 / (N * V * Q - 1)

  # Project subject-specific maps onto design matrix to obtain initial values for
  # S0 and the covariate effects
  Xaug <- cbind(1, X)
  proj = solve(t(Xaug) %*% Xaug) %*% t(Xaug)

  S0   <- matrix(0, nrow = V, ncol = Q)
  Beta <- array(0, dim = c(V, P*Q))

  sigma_2_sq <- rep(0, Q)

  ic_sequence = seq(from = 1, to = N*Q, by = Q)
  cov_sequence = 1:P

  for (iIC in 1:Q){
    est =  Si[,ic_sequence] %*% t(proj)
    S0[, iIC] = est[, 1]
    Beta[, cov_sequence] = est[, 2:(P+1)]

    errs = Si[,ic_sequence] - est %*% t(Xaug)
    sigma_2_sq[iIC] = var(c(errs))

    ic_sequence  <- ic_sequence + 1
    cov_sequence <- cov_sequence + P
  }

  sigma_q_sq = sigma_1_sq + sigma_2_sq

  # DPM Terms
  S0_DPM = list()
  S0_DPM$miu_h      <- rep(0, 150)
  S0_DPM$sigma_sq_h <- rep(0.1, 150)
  S0_DPM$sticks     <- rep(0, 150)
  S0_DPM$n_in_cluster <- rep(0, 150)

  # Split into n_initial_cluster clusters
  # will likely not converge, which is fine. Wrapped to ignore warning message
  suppressWarnings(
    S_kmeans <- kmeans(c(S0), n_initial_cluster, nstart = 5, algorithm="MacQueen")
  )

  # TODO - consider sorting by size?

  # Get DPM information from the clustering results
  S0_DPM$miu_h[1:n_initial_cluster]        <- S_kmeans$centers
  S0_DPM$sigma_sq_h[1:n_initial_cluster]   <- S_kmeans$withinss / (S_kmeans$size)
  S0_DPM$cluster_membership                <- matrix(S_kmeans$cluster, nrow = V, ncol = Q)
  S0_DPM$n_in_cluster[1:n_initial_cluster] <- S_kmeans$size
  S0_DPM$sticks[1:n_initial_cluster]       <- cumprod(runif(10))
  S0_DPM$u                                 <- matrix(runif(Q*V,
                                                           min = 0.0,
                                                           max = S0_DPM$sticks[S0_DPM$cluster_membership]),
                                                     nrow = V,
                                                     ncol = Q
  )


  # Initialize parameters for Horseshoe prior
  Horseshoe = list()
  Horseshoe$tau_sq    <- matrix(runif(Q*P), nrow = Q, ncol = P)
  Horseshoe$xi        <- matrix(runif(Q*P), nrow = Q, ncol = P)
  Horseshoe$lambda_sq <- array(runif(Q*V*P), dim = c(V, P*Q))
  Horseshoe$nu        <- array(runif(Q*V*P), dim = c(V, P*Q))

  guess = list(
    S0 = S0,
    Si = Si,
    beta = Beta,
    A = A,
    S0_DPM = S0_DPM,
    Horseshoe = Horseshoe,
    sigma_sq_q = sigma_q_sq
  )

  return(guess)

}




