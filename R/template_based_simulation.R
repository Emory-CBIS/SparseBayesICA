# Some functions related to simulation neuroimaging data based on ICA template images

generate_mixing_matrix <- function(Q, gentype = "Normal"){

  # Populate with random values as a starting point
  mix_mat_unnormalized <- matrix(rnorm(Q^2), nrow=Q)

  # Normalize to orthonormal columns
  mix_mat_normalized <-  t(mix_mat_unnormalized) %*%
    sqrtm(solve(mix_mat_unnormalized %*% t(mix_mat_unnormalized)))

  return(mix_mat_normalized)

}


# cov types - 0 or 1 for continuous, 2... for categorical with number of levels
generate_model_matrix <- function(N, cov_types, interactions){



}
