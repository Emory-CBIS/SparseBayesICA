#' Load a specified brain mask
#'
#' @param mask_path The file path to the mask file, which should be in nifti format.
#' @return An R list containing the indices of the voxels in the brain, the total number of voxels,
#' and the mask nifti object.
#' @examples
#'
#' load_mask("Documents/masks/MNI_2mm_mask.nii")
load_mask <- function(mask_path){

  # Check input
  if (!file.exists(mask_path)){
    fileName <- basename(mask_path)
    dirName <- dirname(mask_path)
    stop(paste0("Specified mask file, ", fileName, ", not found in directory ", dirName))
  }

  # Load the mask
  mask <- readNifti(mask_path)

  # Get voxels included in mask
  valid_voxels <- which(mask != 0.0)

  V = length(valid_voxels)

  return(list(mask = mask,
              valid_voxels = valid_voxels,
              V = V))

}



