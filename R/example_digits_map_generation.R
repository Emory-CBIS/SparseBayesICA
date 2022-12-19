#' Internal function to generate a random Q x Q orthonormal matrix for a set of N subjects
#' @keywords internal
generate_mixing_matrix <- function(N, Q){
  A = array(0, dim = c(Q, Q, N))
  for (i in 1:N){
    # random temporary mixing matrix
    Ai <- matrix(runif(Q^2, -3, 3), nrow = Q)
    # make orthonormal
    Ai <- Ai %*% solve(sqrtm(t(Ai) %*% Ai))
    # store
    A[,,i] <- Ai
  }
  return(A)
}

#' Internal function to generate a set of covariate values
#' @keywords internal
generate_covariate_values <- function(N, n_discrete, n_continuous, interactions, coding="reference"){
  n_interaction <- nrow(interactions)
  X <- matrix(0, nrow = N, ncol = n_discrete + n_continuous + n_interaction)

  if (n_discrete > 0){
    for (p in 1:n_discrete){
      # make sure not all values for the covariate are the same
      # going to use restriction that at least 30% must be different
      valid_draw <- FALSE
      while(!valid_draw){
        xp <- rbinom(N, 1, prob = 0.5)
        ncat1 <- sum(xp == 1)
        if ( (ncat1 < ceiling(N/3)) | ( (N - ncat1) < ceiling(N/3)) ){
          valid_draw = FALSE
        } else {
          valid_draw = TRUE
        }
      }
      # Check if effects coding requested
      if (coding == "effects"){
        xp <- 2*xp - 1
      }
      X[,p] <- xp
    }
  }

  if (n_continuous > 0){
    # Draw the continous covariate settings
    for (p in (1+n_discrete):(n_discrete + n_continuous)){
      xp <- rnorm(N)
      xp <- (xp - mean(xp)) / sd(xp)
      X[,p] <- xp
    }
  }

  # Create the interaction terms
  X_main_effects <- X[,1:(n_discrete + n_continuous)]
  if (n_interaction > 0){
    for (p in 1:n_interaction){
      xp <- apply(X_main_effects[,which(interactions[p, ] == 1)], 1, prod)
      X[,p+n_discrete+n_continuous] <- xp
    }
  }
  return(X)
}


#' Function to draw a subject specific maps for the digits example
#' @keywords internal
generate_Si <- function(S0, beta, X, sigma_sq_q){
  N <- nrow(X)
  P <- ncol(X)
  Q = nrow(S0)
  V = ncol(S0)
  Si_all <- array(0, dim=c(Q, V, N))
  for (i in 1:N){
    Si <- S0
    for (p in 1:P){
      Si <- Si + X[i, p] * beta[,,p]
    }
    # Add the noise
    for (q in 1:Q){
      Si[q, ] = Si[q, ] + rnorm(V, sd = sqrt(sigma_sq_q[q]))
    }
    Si_all[,,i] <- Si
  }
  return(Si_all)
}



#' Function to create fMRI time courses for the digits example
#' @keywords internal
generate_time_courses <- function(A, Si, sigma_sq_y){
  N <- dim(A)[3]
  Q <- dim(Si)[1]
  V <- dim(Si)[2]
  Yi_all <- array(0, dim = c(Q, V, N))
  for (i in 1:N){
    Yi <- A[,,i] %*% Si[,,i] + matrix(rnorm(Q*V, sd = sqrt(sigma_sq_y)), nrow=Q)
    # Tweak so that generative model is true
    YYt <- Yi %*% t(Yi)
    ed <- eigen(YYt)
    A[,,i] <- t(ed$vectors) %*% A[,,i]
    Yi <- t(ed$vectors) %*% Yi
    Yi_all[,,i] <- Yi
  }
  return(list(A = A,
              Y = Yi_all))
}


#' Function to create fMRI time courses for the digits example
#' @keywords internal
write_data_to_TxV <- function(Y, nT, output_directory, mask_info, use_Di = FALSE){

  Q = dim(Y)[1]
  V = dim(Y)[2]
  N = dim(Y)[3]

  niifiles <- NULL

  tsfm_all       <- array(0, dim = c(Q, nT, N))
  covariance_all <- array(0, dim = c(nT, nT, N))

  for (i in 1:N){
    # Generate a random covariance matrix
    covariance <- LaplacesDemon::rinvwishart(nT + 1, diag(rep(1, nT)))
    covariance_all[,,i] <- covariance

    # Get the eigenvectors that give rotation to this covariance structure
    ed <- eigen(covariance)
    E <- ed$vectors[,1:Q]

    tsfm <- E

    # Rotate the data and increase the dimension
    Yi_TxV <- tsfm %*% rbind(Y[,,i])
    Yi_TxV <- Yi_TxV + matrix(rnorm(nT*V, sd = 0.001), nrow = nT, ncol = V)

    tsfm_all[,,i] <- tsfm

    # convert to a Nifti file with T time points
    Yi_nii <- array_to_nifti(t(Yi_TxV), mask_info, reference_image = NULL)

    # write file to data directory
    fname <- file.path(output_directory, paste0("subject_", i, ".nii"))
    writeNifti(Yi_nii, file = fname)
    niifiles <- c(niifiles, fname)
  }
  return( list(niifiles = niifiles,
               tsfm = tsfm_all,
               trueCov = covariance_all) )
}




#' Function to draw a spatial map with digits as the activation region
#' @keywords internal
draw_shape <- function(multvec, xoff, yoff, block_size, block_thick, imgsize){

  img <- matrix(0, nrow = imgsize, ncol = imgsize)

  # Vertical pieces
  # lower left
  img[(1):(block_thick) + xoff,
      (block_thick + 1 - floor(block_thick)/2 ):(block_thick + block_size + floor(block_thick)/2) + yoff] <- multvec[1]
  # upper left
  img[(1):(block_thick) + xoff,
      (block_thick + 1  - floor(block_thick)/2 ):(block_thick + block_size + floor(block_thick)/2) + yoff + block_size + block_thick] <- multvec[2]
  # lower right
  img[ block_thick + (block_size + 1):(block_size + block_thick) + xoff,
       (block_thick + 1 - floor(block_thick)/2 ):(block_thick + block_size + floor(block_thick)/2) + yoff] <- multvec[3]
  # upper right
  img[ block_thick + (block_size + 1):(block_size + block_thick) + xoff,
       (block_thick + 1 - floor(block_thick)/2 ):(block_thick + block_size + floor(block_thick)/2) + yoff + block_size + block_thick] <- multvec[4]

  # 3 Horizontal pieces
  # bottom
  img[(1+block_thick):(block_thick+block_size) + xoff,
      (1):(block_thick) + yoff] <- multvec[5]
  # middle
  img[(1+block_thick):(block_thick+block_size) + xoff,
      (1):(block_thick) + yoff + block_size + block_thick] <- multvec[6]
  # top
  img[(1+block_thick):(block_thick+block_size) + xoff,
      (1):(block_thick) + yoff + 2*block_size + 2*block_thick] <- multvec[7]

  return(img)
}

#' Function to draw a spatial map with numbers as the activation region
#' @keywords internal
get_shape <- function(number, xoff, yoff, block_size = 10, block_thick = 5, imgsize = 100){


  multvec <- rep(1, 7)

  switch(number,
         "0" = {multvec[6] <- 0},
         "1" = {multvec[1:2] <- 0; multvec[5:7] <- 0},
         "2" = {multvec[2] <- 0; multvec[3] <- 0},
         "3" = {multvec[1:2] <- 0},
         "4" = {multvec[c(1,  5, 7)] <- 0},
         "5" = {multvec[c(1, 4)] <- 0}
  )

  img <- draw_shape(multvec, xoff, yoff, block_size, block_thick, imgsize)

  return(
    list(
      multvec = multvec,
      img = img
    )
  )

}

#' Generate Digits Example Data
#'
#' @details
#' This function will generate a dataset for demonstrating the SparseBayes ICA
#' approach. The data consist of a single slice, with activation regions in the
#' shape of digits from 0-9. Covariate effects are placed within these activation
#' regions.
#'
#' @param N the number of subjects
#' @param Q the number of independent components
#' @param sigma_sq_y scalar variance of the time series noise
#' @param sigma_sq_q vector of length Q containing the between subject variances in the components
#' @param population_map_mean scalar average for the population level map activation regions
#' @param population_map_var scalar variance for the population level map activation regions
#' @param beta_mean scalar average for the population level map activation regions
#' @param beta_var scalar variance for the population level map activation regions
#' @param n_continuous_cov number of continuous covariates in the data. Default is 1
#' @param n_discrete_cov number of categorical covariates in the data (e.g. treatment group). Default is 0.
#' @param interactions matrix of interactions. Each row corresponds to a single
#' interaction. The columns are 1 if the variable is included in that interaction
#' and 0 otherwise. Default is NA, which means no interactions in the model.
#' @param slice_width The number of voxels along one side of the slice. The total
#' number of voxels in the brain will be slice_width squared.
#' @param time_points the number of time points to generate for each subject

#' @return
#' \itemize{
#'   \item data_directory - where the data were stored
#'   \item A - the mixing matrix for each subject
#'   \item S0 - the population level maps
#'   \item Si - the subject level maps
#'   \item beta - the covariate effect spatial maps
#'   \item Y - Q x V time series data for each subject
#'   \item X - matrix with the covariates for each subject
#'   \item sigma_sq_q - same as argument
#'   \item sigma_sq_y - same as argument
#' }
#'
#' @export
generate_digits_example_data <- function(N,
                                          Q,
                                         sigma_sq_y,
                                         sigma_sq_q,
                                         population_map_mean,
                                         population_map_var,
                                         beta_mean,
                                         beta_var,
                                          n_continuous_cov = 1,
                                          n_discrete_cov = 0,
                                          interactions = NA,
                                          slice_width = 75,
                                         time_points = 100){

  # Make sure user input something this function can handle
  if (Q > 10){
    stop("This data generation function can only draw maps for 10 or fewer ICs.")
  }
  if (Q < 3){
    stop("At least 3 components must be used.")
  }
  if (N < 5){
    stop("Please use a sample size of at least 5.")
  }

  #

  # Get the output path (if not specified)
  output_directory <- tools::R_user_dir("SparseBayesICA", which = "data")
  output_directory <- file.path(output_directory, "digits")

  # Create the output directory
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

  # Determine the number of voxels using the slice_width
  V <- slice_width^2

  # Determine the total number of covariates
  if (is.na(interactions)){
    interactions <- array(dim = c(0,n_continuous_cov + n_discrete_cov))
  }
  P <- n_continuous_cov + n_discrete_cov + nrow(interactions)

  # Some calculations to enable drawing the numbers
  block_size = 10
  block_thick = 5
  maximum_starting_index <- slice_width - 3.5 * block_size

  # Generate the number for each IC
  S0   <- matrix(0, nrow = Q, ncol  = V)
  beta <- array(0, dim = c(Q, V, P))
  for (q in 1:Q){

      xoffset <- sample(1:maximum_starting_index, size = 1)
      yoffset <- sample(1:maximum_starting_index, size = 1)

      shape_info <- get_shape(paste(q-1), xoffset, yoffset, imgsize = slice_width)

      # S0
      S0[q,] <- as.vector(shape_info$img) * (population_map_mean + sqrt(population_map_var) * rnorm(V))

      # Betas (overlap with shape)
      for (p in 1:P){
        # Pick a random part of the shape to edit with the covariate effect
        indselect       <- sample(which(shape_info$multvec == 1), 1)
        # This controls which parts of the shape are drawn
        beta_multvec    <- rep(0, 7)
        beta_multvec[indselect] <- 1
        # Draw the shape
        beta_shape_info <- draw_shape(beta_multvec,
                                      xoffset,
                                      yoffset,
                                      block_size = block_size,
                                      block_thick = block_thick,
                                      imgsize = slice_width)
        # Store the shape
        beta[q,,p] <- as.vector(beta_shape_info) * (beta_mean + sqrt(beta_var) * rnorm(V))
      }
  }

  # Create a "mask"
  mask_info <- list()
  mask_info$valid_voxels <- 1:V
  mask_info$mask <- array(1, dim =c(slice_width, slice_width, 1))

  # Save the mask to the output directory
  mask_final <- matrix_to_nifti(c(mask_info$mask), mask_info, reference_image = NULL)
  fname <- file.path(output_directory, "mask.nii")
  writeNifti(mask_final, file = fname)

  # Generate the mixing matrices for each subject
  A_initial <- generate_mixing_matrix(N, Q)

  # Get the covariates for each subject
  X <- generate_covariate_values(N, n_discrete_cov, n_continuous_cov, interactions)

  # Generate the subject-specific maps
  Si <- generate_Si(S0, beta, X, sigma_sq_q)

  # Generate the preprocessed time course data
  tc_result <- generate_time_courses(A_initial, Si, sigma_sq_y)
  Y <- tc_result$Y
  A <- tc_result$A

  # Color the data to T x V and write out nifti files
  tsfm_list <- write_data_to_TxV(Y, time_points, output_directory, mask_info, use_Di = use_Di)

  niifiles = tsfm_list$niifiles
  tsfm     = tsfm_list$tsfm
  trueCov  = tsfm_list$trueCov

  # Write out the covariates with the nifti filenames
  X_covcol <- cbind(niifiles, X)
  fname <- file.path(output_directory, "covariates.csv")
  write.csv(X_covcol, file = fname, row.names = FALSE)

  # Save .Rdata file with all simulation quantities
  fname <- file.path(output_directory, "truth.RData")
  true_values <- list(
    data_directory = output_directory,
    A = A,
    S0 = S0,
    Si = Si,
    beta = beta,
    Y = Y,
    X = X,
    sigma_sq_q = sigma_sq_q,
    sigma_sq_y = sigma_sq_y
  )
  save(true_values,
       file = fname)

  return(true_values)
}

