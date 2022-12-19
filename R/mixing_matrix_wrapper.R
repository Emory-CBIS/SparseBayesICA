#' wrapper for mixing matrix function
#' @keywords internal
sample_mixing_matrix <- function(A,
                                 YSi,
                                 YSiOmega,
                                 Omega,
                                 YYt,
                                 angle_range,
                                 sin_theta,
                                 cos_theta){


  # Number of subjects (potentially cluster specific if in parallel)
  Q <- dim(A)[2]
  N <- dim(A)[1] / Q

  endIndex = 0
  negYSiSigInvover2 <- -Omega/2

  for (i in 1:N){
    startIndex = endIndex + 1
    endIndex   = endIndex + Q

    A[startIndex:endIndex,] <- .gibbs_sample_mixing_matrix(A[startIndex:endIndex,],
                                                         (YSiOmega[startIndex:endIndex,]),
                                                         (negYSiSigInvover2),
                                                         (YYt[startIndex:endIndex,]),
                                                         angle_range, sin_theta,
                                                         cos_theta,
                               column_reduced_mixing_matrix,
                               null_space)


    if (any(is.na(A[startIndex:endIndex,]))){
      print(paste0("Subject ", i))
      stop("NAs in mixing matrix")
    }

  }

  return(A)
}





#' wrapper for mixing matrix function
#' @keywords internal
sample_mixing_matrix_inplace <- function(A,
                                 YSi,
                                 YSiOmega,
                                 Omega,
                                 YYt,
                                 angle_range,
                                 sin_theta,
                                 cos_theta){


  # Number of subjects (potentially cluster specific if in parallel)
  Q <- dim(A)[2]
  N <- dim(A)[1] / Q

  endIndex = 0
  negYSiSigInvover2 <- -Omega/2

  for (i in 1:N){
    startIndex = endIndex + 1
    endIndex   = endIndex + Q

    gibbs_sample_mixing_matrix_inplace(A[startIndex:endIndex,],
                                                          (YSiOmega[startIndex:endIndex,]),
                                                          (negYSiSigInvover2),
                                                          (YYt[startIndex:endIndex,]),
                                                          angle_range, sin_theta,
                                                          cos_theta,
                                                          column_reduced_mixing_matrix,
                                                          null_space)


    if (any(is.na(A[startIndex:endIndex,]))){
      print(paste0("Subject ", i))
      stop("NAs in mixing matrix")
    }

  }

}

