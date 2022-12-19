#' Load a brain mask from a Nifti file
#'
#' @param mask_path The file path to the mask file, which should be in nifti format.
#' @return A list containing the following:
#' \itemize{
#'   \item mask - the nifti object for the mask
#'   \item valid_voxels - a list containing the linear indices of the voxels that
#'   are in the masked region
#'   \item V - the total number of voxels in the brain mask
#' }
#'
#' @export
load_mask <- function(mask_path){

  # Check input
  if (!file.exists(mask_path)){
    fileName <- basename(mask_path)
    dirName <- dirname(mask_path)
    stop(paste0("Specified mask file, ", fileName, ", not found in directory ", dirName))
  }

  # Load the mask
  mask <- readNifti(mask_path)

  # Fix potential dimension problem for 2d slices
  if (length(dim(mask)) == 2){
    dim(mask) <- c(dim(mask)[1], dim(mask)[2], 1)
  }

  # Get voxels included in mask
  valid_voxels <- which(mask != 0.0)

  V = length(valid_voxels)

  return(list(mask = mask,
              valid_voxels = valid_voxels,
              V = V))

}



