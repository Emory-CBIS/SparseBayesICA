## code to prepare `example` dataset goes here


# TEMPORARY
td <- "/Users/joshlukemire/Documents/Research/cbis/SparseBayesICA/revision_simulation/templates"
mask <- load_mask(file.path(td, "lowres_mask.nii"))

example_RSNs   <- matrix(0, nrow = 21718, ncol = 10)
example_covariate_effects  <- matrix(0, nrow = 21718, ncol = 10)
for (ic in 1:10){
  x <- readNifti(file.path(td, paste0("rsn_", ic, ".nii")))
  vv <-x[mask$valid_voxels]
  example_RSNs[,ic] <- vv

  x <- readNifti(file.path(td, paste0("rsn_", ic, "_beta_1_small.nii")))
  vv <-x[mask$valid_voxels]
  example_covariate_effects[,ic] <- vv
}

usethis::use_data(example_RSNs, example_covariate_effects,
                  internal = TRUE,
                  overwrite = TRUE)
