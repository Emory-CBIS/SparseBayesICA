#' Prewhiten a subject's fMRI time courses
#'
#' @details
#' `preprocess_subject` takes a subject's fMRI time courses and applies a
#' prewhitening procedure.
#'
#' @param data A T x V dataset where T is the number of fMRI volumes and V is the
#' number of voxels.
#' @param Q The requested dimension of the prewhitened data. This should be equal
#' to the number of independent components that you plan to include in the model.
#' @return A list containing the prewhitnened data and the eigenvalues of the
#' original scale data covariance
#'
#' @examples
#' library(SparseBayesICA)
#'
#' # Get the filepath to the example data mask (provided in this package)
#' mask_path  <- system.file("extdata", "example_mask.nii",
#'                            package = "SparseBayesICA", mustWork = TRUE)
#'
#' # Load the mask
#' mask_info  <- load_mask(mask_path)
#'
#' # Filepath for one of the example subjects provided in this package
#' nifti_path <- system.file("extdata", "example_subject_1.nii",
#'                            package = "SparseBayesICA", mustWork = TRUE)
#'
#' # Load the data
#' nifti_data <- readNifti(nifti_path)
#'
#' # Apply the mask to the data to get a T x V data matrix instead of 4D data
#' masked_data <- nifti_to_matrix(nifti_data, mask_info$valid_voxels)
#'
#' # Prewhiten this subject's data to Q = 5 whitened time points
#' prewhitened_results <- preprocess_subject(masked_data, Q = 5)
#'
#' @seealso [load_mask()]
#' @export
preprocess_subject <- function(data, Q){

  # Dimensions
  Ti = dim(data)[1]
  V  = dim(data)[2]

  # Verify input
  if (Ti < Q){
    stop(paste0("Dimension reduction to ", Q, " points requested but subject's",
                " data only contains ", Ti, " time points."))
  }
  if (V < Ti){
    warning("Number of voxels is less than number of time points Check if data matrix",
         " needs to be transposed.")
  }

  # Center voxels, this is done in case previous processing gave them some other mean
  # (for example, mean = 100 is common output from some pipelines)
  data <- sweep(data, 2, colMeans(data))

  # Center time points
  data <- sweep(data, 1, rowMeans(data))

  # PCA dimension reduction
  data_pca <- prcomp(t(data), center = FALSE)

  # Get the eigenvals
  eigvals_Q <- data_pca$sdev[1:Q]^2
  mean_remaining_eigvals <- mean(data_pca$sdev[-c(1:Q)]^2)

  # Eigenvectors
  U_eigen <- data_pca$rotation[, 1:Q]

  # Construct the transformation matrix
  tsfm_matrix <- t(U_eigen)
  tsfm_matrix <- diag( (eigvals_Q)^(-1/2) ) %*% tsfm_matrix

  # Construct the whitened data and flip it to be V x Q
  Y = t(tsfm_matrix %*% data)

  # Make sure the whitened data are centered (time courses mean 0)
  Y <- sweep(Y, 2, colMeans(Y))

  prewhitening_results <- list(Y = Y, eigenvalues = data_pca$sdev^2)

  # Returns whitened data
  return(prewhitening_results)
}

#' Dimensionality reduction of a set of fMRI time courses to a set of principal components
#'
#' @description
#' `PCA_dimension_reduction()` takes a set of fMRI time courses and reduces it
#' to a set of `nPC` orthogonal time courses. This is useful for creating
#' dimension reduced subject and group level data for temporal concatenation
#' group ICA, which can be used to obtain starting values for SparseBayes ICA
#' if desired.
#'
#' @details
#' Let \eqn{Y} be the fMRI time courses with T rows and V columns. After centering,
#' we can apply a spectral decomposition to the \eqn{T \times T} data covariance
#' to obtain:
#' \deqn{Cov(Y) = EDE'}
#' This function then returns \eqn{E^* D^*}, where \eqn{E^*} consists of the
#' first `nPC` columns of \eqn{E} and \eqn{D^*} is a diagonal matrix with the
#' `nPC` largest eigenvalues of \deqn{Cov(Y)}.
#'
#' @param data A T x V dataset where T is the number of fMRI volumes and V is the
#' number of voxels.
#' @param nPC The requested number of principal components
#' @return A nPC x V data matrix with rows containing the principal components
#' @examples
#' library(SparseBayesICA)
#'
#' # Get the filepath to the example data mask (provided in this package)
#' mask_path  <- system.file("extdata", "example_mask.nii",
#'                            package = "SparseBayesICA", mustWork = TRUE)
#'
#' # Load the mask
#' mask_info  <- load_mask(mask_path)
#'
#' # Filepath for one of the example subjects provided in this package
#' nifti_path <- system.file("extdata", "example_subject_1.nii",
#'                            package = "SparseBayesICA", mustWork = TRUE)
#'
#' # Load the data
#' nifti_data <- readNifti(nifti_path)
#'
#' # Apply the mask to the data to get a T x V data matrix instead of 4D data
#' masked_data <- nifti_to_matrix(nifti_data, mask_info$valid_voxels)
#'
#' # Perform dimension reduction to 10 principal components
#' subject_PCA <- PCA_dimension_reduction(masked_data, nPC = 10)
#'
#' @export
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







#' Use fastica to obtain a set of initial values for the population level components.
#'
#' @description
#' `obtain_s0_guess_fastica()` takes a set of subject level time courses that
#' have been stacked across the time domain and reduced to dimension Q. It then
#' uses fastica to obtain the group level components, which can be used as an
#' initial guess for the population level components in the SparseBayes ICA model.
#'
#' @param data A Q x V dataset that was obtained by performing PCA on the stacked,
#' dimension reduced subject level data. See \code{vignette("digits_example", package = "SparseBayesICA")}
#' for more details.
#' @return A V x Q matrix, where each column contains a starting value
#' for an independent component.
#' @export
obtain_s0_guess_fastica <- function(data){

  Q = dim(data)[1]
  V = dim(data)[2]

  # Apply fastica
  ica_decomposition <- fastICA(t(data), Q)

  S0 = ica_decomposition$S

  return(S0)

}








#' Use the population level components to obtain a set of initial values for all
#' SparseBayes ICA model parameters
#'
#' @description
#' `obtain_initial_guess_sparsebayes()` Uses the population level initial values
#' generated from `obtain_s0_guess_fastica()` to obtain initial parameter settings
#' for all SparseBayes ICA model parameters.
#'
#' @param data A list containing the prewhitened data for each subject.
#' @param S0 The initial guess for the population level parameters
#' @param X A matrix of covariates with N rows and P columns
#' @param n_initial_cluster The number of initial clusters for the DPM. Default is 10.
#' @return A list object containing the initial parameters for SparseBayes ICA
#' @export
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


#' Get an initial guess using a random set of starting values
#'
#' @description
#' `random_initialization()` Generates a random set of starting values for
#' SparseBayes ICA. In general, this should be avoided in favor of using TCGICA
#' and dual regression.
#'
#' @param data A list containing the prewhitened data for each subject.
#' @param X A matrix of covariates with N rows and P columns
#' @param n_initial_cluster The number of initial clusters for the DPM. Default is 10.
#' @return A list object containing the initial parameters for SparseBayes ICA
random_initialization <- function (data, X, n_initial_cluster = 10){

  # Dimensions
  N = nrow(X)
  V = dim(data[[1]])[1]
  Q = dim(data[[1]])[2]
  P = ncol(X)

  S0 <- matrix(runif(Q*V, min = -1, max = 1), nrow = V)

  Si = array(0, dim=c(V, N*Q) )
  A  = array(0, dim=c(N*Q, Q) )

  sse_level1 <- 0.0
  eI = 0
  for (i in 1:N){

    # Project
    Atemp = matrix(rnorm(Q^2), nrow = Q)

    # Apply orthonormal tsfm
    sI = eI + 1
    eI = eI + Q
    A[sI:eI,] = Atemp %*% sqrtm(solve(t(Atemp)%*%Atemp));

    # Back reconstruct Si using new A
    Si[,sI:eI] = matrix(rnorm(Q*V), nrow = V)

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
    S0[, iIC] = runif(V, min = -1, max = 1)
    Beta[, cov_sequence] = runif(V*P, min = -1, max = 1)

    sigma_2_sq[iIC] = 1.0

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




