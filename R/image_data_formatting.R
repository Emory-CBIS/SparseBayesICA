# Functions related to NIFTI to array manipulation

nifti_to_matrix <- function(data, valid_voxels){

  V = length(valid_voxels)

  if (length(dim(data)) == 3){
    dim(data)[4] <- 1
  }
  Ti = dim(data)[4]

  matrix_data <- matrix(0, nrow = Ti, ncol = V)

  for (t in 1:Ti){
    data_t = data[,,,t]
    matrix_data[t, ] = data_t[valid_voxels]
  }

  return(matrix_data)

}

matrix_to_nifti <- function(data, mask_info, reference_image = NULL){

  new_image <- array(0, dim = dim(mask_info$mask))

  new_image[mask_info$valid_voxels] <- data

  if (is.null(reference_image)){
    reference_image = mask_info$mask
  }

  new_nifti <- asNifti(new_image, reference = reference_image)

  return(new_nifti)

}

# Third dimension should be volumes or time points
array_to_nifti <- function(data, mask_info, reference_image = NULL){

  # TODO check 3d, otherwise suggest matrix_to_nifti function

  nVol = dim(data)[2]

  temp_image <- array(0, dim = c(dim(mask_info$mask)))
  new_image <- array(0, dim = c(dim(mask_info$mask), nVol))

  for (iVol in 1:nVol){
    temp_image[mask_info$valid_voxels] <-  t(data[,iVol])
    new_image[,,,iVol] <- temp_image
  }

  if (is.null(reference_image)){
    reference_image = mask_info$mask
  }

  new_nifti <- asNifti(new_image, reference = reference_image)

  return(new_nifti)

}


# Function to create a tibble from the sparse bayes findings that can be used to
# plot images
flatten_sparsebayes_results <- function(results, mask_info, covariate_labels = NULL){

  # Create the basic part of the frame
  base_layout <- reshape2::melt(mask_info$mask)
  base_layout <- base_layout[mask_info$valid_voxels,] %>%
    dplyr::select(-c("value"))
  colnames(base_layout) <- c("sag", "cor", "axi")

  # Get dimensions from sparsebayes results
  Q <- ncol(results$S0_posterior_mean)
  P <- dim(results$Beta_posterior_mean)[3]

  if (is.null(covariate_labels)){
    covariate_labels <- paste0("Cov", 1:P)
  } else {
    if (length(covariate_labels) != P){
      stop( paste0("Length of covariate labels, ", length(covariate_labels),
                   " does not match number of covariates, ", P) )
    }
  }

  all_results <- tibble()


  for (q in 1:Q){

    new_rows <- cbind(base_layout, Q = q, S0 = results$S0_posterior_mean[,q])
    new_rows <- cbind(new_rows, results$Beta_posterior_mean[,q,])
    colnames(new_rows)[6:(6+P-1)] <- covariate_labels
    new_rows <- cbind(new_rows, results$Beta_pvals[,q,])
    colnames(new_rows)[(6+P):(6+2*P-1)] <- paste0(covariate_labels, "_pvalue")

    all_results <- rbind(all_results, new_rows)

  }


  return(all_results)

}




