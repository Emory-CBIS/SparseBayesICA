#' Flatten a Nifti image to a matrix
#'
#' @details
#' `nifti_to_matrix` takes a subject's fMRI volumes from a NIFTI file and converts them
#' to a T x V data matrix, where T is the number of volumes and V is the number of
#' voxels in the brain.
#'
#' @param data 4 dimensional Nifti data for a subject
#' @param valid_voxels A list of the linear indices of the voxels in the brain mask.
#' This can be obtained using \code{\link{load_mask}}
#' @return A T x V data matrix
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
#' @seealso [load_mask()]
#' @export
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

#' Take a matrix and use it to create a brain volume
#'
#' @details
#' `matrix_to_nifti` takes a matrix with V (number of voxels) elements and places
#' the values into a brain. The mappings from the matrix to the brain are provided by
#' the mask_info object, which can be obtained using \code{\link{load_mask}}.
#'
#' @param data A data matrix with V elements
#' @param mask_info The output from using \code{\link{load_mask}} on a nifti file
#' containing a brain mask for a brain with V voxels.
#' @param reference_image Imaging containing other header information. Default is NULL, which
#' results in the mask file header being used.
#' @return A nifti volume
#'
#' @seealso [nifti_to_matrix()]
#' @export
matrix_to_nifti <- function(data, mask_info, reference_image = NULL){

  new_image <- array(0, dim = dim(mask_info$mask))

  new_image[mask_info$valid_voxels] <- data

  if (is.null(reference_image)){
    reference_image = mask_info$mask
  }

  new_nifti <- asNifti(new_image, reference = reference_image)

  return(new_nifti)

}


#' Take an array and use it to create a nifti file with multiple volumes
#'
#' @details
#' `array_to_nifti` takes a matrix with V (number of voxels) elements and L columns
#' and places
#' the values into a brain with L volumes. The mappings from the matrix to the brain are provided by
#' the mask_info object, which can be obtained using \code{\link{load_mask}}. This
#' will be renamed eventually to matrix_to_nifti, and that function will be discontinued.
#'
#' @param data A data matrix with V elements and L columns (e.g. ICs)
#' @param mask_info The output from using \code{\link{load_mask}} on a nifti file
#' containing a brain mask for a brain with V voxels.
#' @param reference_image Imaging containing other header information. Default is NULL, which
#' results in the mask file header being used.
#' @return A nifti object with L volumes
#'
#' @seealso [nifti_to_matrix()]
#' @export
array_to_nifti <- function(data, mask_info, reference_image = NULL){

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

