if (0 == 1){


  library(hcicaR)

  # data path
  dp <- "/Users/joshlukemire/Documents/Research/cbis/SparseBayesICA/simulation_demo_data"

  # subject list
  subject_list <- file.path(dp, paste0("subject_", 1:50, ".nii"))

  # mask file
  mask_file <- file.path(dp, "lowres_mask.nii")

  # covariate file
  cov_file <- file.path(dp, "covariates.csv")

  # Covariate data
  X <- read.csv(cov_file)
  N <- nrow(X)
  X <- X[,2:4]
  ff <- rep(1, N) ~ 0 + x1 + x2 + x3
  utils::str(m <- model.frame(ff, X))
  mat <- model.matrix(ff, m)
  X = mat

  # Preprocess subjects
  Q = 5
  nPC = 10

  mask_info = load_mask(mask_file)

  whitened_data <- list()

  stacked_pca_data <- matrix(0, nrow = nPC * N, ncol = mask_info$V)

  for (i in 1:N){
    print(i)
    data = readNifti(file = subject_list[i])
    data_mat <- nifti_to_matrix(data, mask_info$valid_voxels)
    whitened_data[[i]] <-  preprocess_subject(data_mat, Q)

    # Get the contribution to the 2 stage PCA reduction
    ind_start <- (i-1)*nPC + 1
    ind_end <- i*nPC
    stacked_pca_data[ind_start:ind_end, ] <- PCA_dimension_reduction(data_mat, nPC)
  }

  pca_reduced_group_data <- PCA_dimension_reduction(stacked_pca_data, Q)

  # Fast ICA to get initial components
  population_average_maps <- obtain_s0_guess_fastica(pca_reduced_group_data)

  # Initial values for SparseBayes ICA based on population average maps above
  initial_values = obtain_initial_guess_sparsebayes(whitened_data, population_average_maps, X)


  # Define stuff so can test sparse bayes function
  Y = array(0, dim = c(mask_info$V, N*Q))
  qSeq <- 1:Q
  for (i in 1:N){
    Y[,qSeq] = whitened_data[[i]]
    qSeq <- qSeq + Q
  }






















  # Write workspace to R data to speed things along
  library(hcicaR)
  fname = "/Users/joshlukemire/Desktop/testRmats.RData"
  #save(Y, initial_values, X, file = "/Users/joshlukemire/Desktop/testRmats.RData")

  dp <- "/Users/joshlukemire/Documents/Research/cbis/SparseBayesICA/simulation_demo_data"
  mask_file <- file.path(dp, "lowres_mask.nii")
  mask_info = load_mask(mask_file)

  load(fname)

  Q = 5
  N = 50

  nburnin = 100
  nmcmc = 100

  initial_values$sigma_sq_q <- initial_values$sigma_q_sq


  # SWITCH TO THIS!!!
  # https://cran.r-project.org/web/packages/parallelly/readme/README.html

  cl = makeCluster(4, outfile="/Users/joshlukemire/Desktop/output.txt")
  results = SparseBayesICAPar(cl, Y, X, initial_values, nburnin = 1500, nmcmc = 1500)
  # results = SparseBayesICAPar(cl, Y, X, initial_values, nburnin = nburnin, nmcmc = nmcmc)
  # results = SparseBayesICAPar(cl, Y, X, initial_values, nburnin = 1000, nmcmc = 1000)
  #results = SparseBayesICAPar(cl, Y, X, initial_values, nburnin = 20, nmcmc = 40)
  stopCluster(cl)

  initial_values$A[1:5, 1:5]
  results$Ai_posterior_mean[1:5, 1:5, 1]

  cor(initial_values$A[1:5, 1:5],
      results$Ai_posterior_mean[1:5, 1:5, 1])

  cor(initial_values$A[11:15, 1:5],
      results$Ai_posterior_mean[1:5, 1:5, 3])




  # SAve nifti file with results
  S0nii <- matrix_to_nifti(results$S0_posterior_mean[,1], mask_info)
  writeNifti(S0nii, "/Users/joshlukemire/Desktop/testIC1.nii")
  S0nii <- matrix_to_nifti(results$S0_posterior_mean[,2], mask_info)
  writeNifti(S0nii, "/Users/joshlukemire/Desktop/testIC2.nii")
  S0nii <- matrix_to_nifti(results$S0_posterior_mean[,3], mask_info)
  writeNifti(S0nii, "/Users/joshlukemire/Desktop/testIC3.nii")

  # SAve nifti file with results
  beta1_ic1nii <- array_to_nifti(results$Beta_posterior_mean[,,1], mask_info)
  writeNifti(beta1_ic1nii, "/Users/joshlukemire/Desktop/testbeta1IC1.nii")
  beta1_ic2nii <- matrix_to_nifti(results$Beta_posterior_mean[,Q+1], mask_info)
  writeNifti(beta1_ic2nii, "/Users/joshlukemire/Desktop/testbeta1IC2.nii")
  beta1_ic3nii <- matrix_to_nifti(results$Beta_posterior_mean[,2*Q+1], mask_info)
  writeNifti(beta1_ic3nii, "/Users/joshlukemire/Desktop/testbeta1IC3.nii")


  At = initial_values$A[1:5, 1:5]
  YSit = t(Y[,1:5]) %*% initial_values$Si[,1:5]

  angle_range = seq(from = 0.0, to = 2*pi, length.out = 100)
  sin_theta <- sin(angle_range)
  cos_theta <- cos(angle_range)

  sample_mixing_matrix(At, YSit, YSit, angle_range, sin_theta, cos_theta)

}
