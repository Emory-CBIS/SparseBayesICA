# combinations can be either the combinations or the total size
partition_into_chunks <- function(n_node, combinations){

  # Get the total number of units to partition
  if (length(combinations) > 1){
    n_units <- length(combinations)
  } else {
    n_units <- combinations
    combinations <- 1:n_units
  }

  combination_sizes <- floor(n_units / n_node)

  # Account for uneven splits
  n_adj <- n_units - combination_sizes * n_node

  splits <- list()
  indEnd <- 0
  for (iNode in 1:(n_node-n_adj)){
    indStart <- indEnd + 1
    indEnd   <- indEnd + combination_sizes
    splits[[iNode]] <- combinations[indStart:indEnd]
  }

  if (n_adj > 0){
    for (iNode in (n_node-n_adj+1):n_node){
      indStart <- indEnd + 1
      indEnd   <- indEnd + combination_sizes + 1
      splits[[iNode]] <- combinations[indStart:indEnd]
    }
  }

  return(splits)

}

# See: https://stackoverflow.com/questions/31068289/clusterexport-to-single-thread-in-r-parallel
export_to_cluster <- function(cluster, var_name, variable){

  clusterCall(cluster, function(d) {
    assign(var_name, d, pos=.GlobalEnv)
    NULL
  }, variable)

  NULL
}

update_cluster_vector_variable <- function(cluster, var_name, variable){

  clusterCall(cluster, function(d) {
    update_vector_in_place( eval(as.symbol(var_name)) , d)
    NULL
  }, variable)

  NULL
}

update_cluster_matrix_variable <- function(cluster, var_name, variable){

  clusterCall(cluster, function(d) {
    update_matrix_in_place( eval(as.symbol(var_name)) , d)
    NULL
  }, variable)

  NULL
}

# this requires that the variable on the cluster already has the variable
# with name given by row_index_list_name
update_cluster_matrix_variable_rowspec <- function(cluster, var_name, variable, row_index_list_name){

  clusterCall(cluster, function(d) {
    update_matrix_in_place_rowspec( eval(as.symbol(var_name)) , d, eval(as.symbol(row_index_list_name)))
    NULL
  }, variable)

  NULL
}


export_worker_state <- function(fname){

  save(S0, Beta, posterior_mean, posterior_variance, draw,
       AtY, X_with_int, XtX, sigma_q_sq,
       tau_sq, lambda_sq,
       mu_h_div_sigma_sq_h, sigma_sq_h, cluster_memberships,
       file = fname)

  NULL

}


initialize_sparsebayes_nodes <- function(cl, voxel_combinations, subject_combinations,
                                         Yin, Xin, guess, angle_range, sin_theta, cos_theta, S0_DPM){

  clusterEvalQ(cl, {
    library(hcicaR)
  })

  n_node <- length(cl)

  Q = ncol(guess$S0)
  N = nrow(Xin)

  # Row indices for each node for the mixing matrix update
  # this is used to divide the data up by subject when updating the
  # mixing matrices
  mixing_matrix_node_row_inds <- list()
  for (ic in 1:n_node){
    inds <- subject_combinations[[ic]]
    mixing_matrix_node_row_inds[[ic]] <- as.integer(inds)
  }

  for (i in 1:n_node){
    voxels   <- voxel_combinations[[i]]
    subjects <- subject_combinations[[i]] # elements of Q*V dimensional variables

    # Preprocessed Time Courses
    export_to_cluster(cl[i], 'Y', Yin[voxels, ])

    # S0
    export_to_cluster(cl[i], 'S0', guess$S0[voxels, ])

    export_to_cluster(cl[i], 'S0_posterior_mean',  matrix(0, nrow = length(voxels), ncol = dim(guess$S0)[2]))

    export_to_cluster(cl[i], 'Si',  guess$Si[voxels, ] + 0)

    export_to_cluster(cl[i], 'AtY', guess$Si[voxels, ] + 0)

    export_to_cluster(cl[i], 'Ei',  guess$Si[voxels, ] + 0)

    export_to_cluster(cl[i], 'VxQStore', matrix(0, nrow=length(voxels), ncol=Q))


    # Beta
    export_to_cluster(cl[i], 'Beta', guess$beta[voxels, ])
    export_to_cluster(cl[i], 'Beta_posterior_mean', matrix(0, nrow = dim(guess$beta[voxels, ])[1], ncol = dim(guess$beta)[2]))
    export_to_cluster(cl[i], 'Beta_ngt0',          matrix(0, nrow = dim(guess$beta[voxels, ])[1], ncol = dim(guess$beta)[2]))

    lambda_sq_posterior_mean = matrix(0, nrow = dim(guess$beta[voxels, ])[1], ncol = dim(guess$beta)[2])
    export_to_cluster(cl[i], 'lambda_sq_posterior_mean', lambda_sq_posterior_mean)


    # Ai
    export_to_cluster(cl[i], 'A', guess$A[subjects, ] + 0)
    export_to_cluster(cl[i], 'mixing_matrix_node_row_inds', mixing_matrix_node_row_inds[[i]])
    export_to_cluster(cl[i], 'YSiOmegaChunk', matrix(0, nrow = length(mixing_matrix_node_row_inds[[i]]), ncol = Q ))
    export_to_cluster(cl[i], 'negYSiSigInvover2', matrix(0, nrow = Q, ncol = Q))

    YYt <- matrix(0, nrow=length(subjects), Q)
    for (iSubj in 1:(length(subjects)/Q)){
      inds <-  (Q*(iSubj-1)+1):(Q*iSubj)
      YYt[ inds, ] = t(Yin[,subjects[inds]]) %*% Yin[,subjects[inds]]
    }
    export_to_cluster(cl[i], 'YYt', YYt)

    # DPM u
    export_to_cluster(cl[i], 'u', S0_DPM$u[voxels, ])
    export_to_cluster(cl[i], 'cluster_memberships', S0_DPM$cluster_membership[voxels, ])

    # Horseshoe
    export_to_cluster(cl[i], 'lambda_sq', guess$Horseshoe$lambda_sq[voxels, ])
    export_to_cluster(cl[i], 'nu',  matrix(0.1, nrow = length(voxels), ncol = ncol(guess$Horseshoe$lambda_sq)) )
    export_to_cluster(cl[i], 'cauchy_mixing', matrix(0.1, nrow = length(voxels), ncol = ncol(guess$Horseshoe$lambda_sq)) )
    export_to_cluster(cl[i], 'cauchy_mixing_prior',  matrix(0.1, nrow = length(voxels), ncol = ncol(guess$Horseshoe$lambda_sq)) )


    # Storage space for mixing matrix update variables
    export_to_cluster(cl[i], 'column_reduced_mixing_matrix', matrix(0, nrow = Q, ncol = Q - 2))
    export_to_cluster(cl[i], 'null_space', matrix(0, nrow=Q, ncol=2))



  }

  # Quantities where the entire item is needed on all nodes
  export_to_cluster(cl, 'X', X)
  export_to_cluster(cl, 'angle_range', angle_range)
  export_to_cluster(cl, 'sin_theta', sin_theta)
  export_to_cluster(cl, 'cos_theta', cos_theta)
  export_to_cluster(cl, 'sigma_sq_q', guess$sigma_sq_q + 0)
  export_to_cluster(cl, 'Q', Q)
  export_to_cluster(cl, 'sticks', S0_DPM$sticks + 0)
  export_to_cluster(cl, 'mu_h', S0_DPM$miu_h + 0)
  export_to_cluster(cl, 'sigma_sq_h', S0_DPM$sigma_sq_h + 0)
  export_to_cluster(cl, 'n_in_cluster', S0_DPM$n_in_cluster + 0)
  export_to_cluster(cl, 'tau_sum', matrix(0, nrow = Q, ncol = ncol(X)))
  export_to_cluster(cl, 'Ai', guess$A + 0)
  export_to_cluster(cl, 'YSi', array(0, dim = c(Q*N, Q)))


  # Model matrix
  X_with_int <- matrix(cbind(1, X), nrow = nrow(X))
  XtX <- t(X_with_int) %*% X_with_int
  export_to_cluster(cl, 'X_with_int', X_with_int)
  export_to_cluster(cl, 'XtX',  matrix(XtX, nrow=ncol(X_with_int)))
  export_to_cluster(cl, 'Ip', diag(1, ncol(X_with_int)))


  # Intermediate quantities used in spatial map update
  draw <- rep(0, ncol(X_with_int))
  posterior_mean  <- rep(0, ncol(X_with_int))
  posterior_variance <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  #mapping <- as.integer(rep(0, 150))
  export_to_cluster(cl, 'draw', draw)
  export_to_cluster(cl, 'posterior_mean', posterior_mean)
  export_to_cluster(cl, 'posterior_variance',  matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int)))
  export_to_cluster(cl, 'prior_precision',     matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int)))


  # Sufficient statistics for the DPM means and variances
  export_to_cluster(cl, 'sum_Ssq_qh', matrix(0, nrow = 150, ncol = Q))
  export_to_cluster(cl, 'Sbar_qh', matrix(0, nrow = 150, ncol = Q))
  export_to_cluster(cl, 'n_in_cluster_qh', matrix(0, nrow = 150, ncol = Q))

  # Variance update
  export_to_cluster(cl, 'sse', rep(0, Q))
  export_to_cluster(cl, 'reg_scale', rep(0, Q))

  # Horseshoe
  export_to_cluster(cl, 'tau_sq', guess$Horseshoe$tau_sq + 0)

  # Storage space for mixing matrix update variables
  #export_to_cluster(cl, 'column_reduced_mixing_matrix?', matrix(0, nrow = Q, ncol = Q - 2))
  # export_to_cluster(cl, 'selected_columns_mixing_matrix', matrix(0, nrow = Q, ncol = 2))
  # export_to_cluster(cl, 'selected_columns_proj_null_space', matrix(0, nrow = 2, ncol = 2))
  # export_to_cluster(cl, 'btilde', matrix(0, nrow = 2, ncol = 2))
  # export_to_cluster(cl, 'atilde', matrix(0, nrow = 2, ncol = 2))
  # export_to_cluster(cl, 'rotation_matrix', matrix(0, nrow = 2, ncol = 2))
  # export_to_cluster(cl, 'null_rotation', matrix(0, nrow = Q, ncol = 2))
  # export_to_cluster(cl, 'log_probs', rep(0, 2*length(angle_range)))
  # export_to_cluster(cl, 'probs', rep(0, 2*length(angle_range)))
  #export_to_cluster(cl, 'null_space', matrix(0, nrow=Q, ncol=2))

  export_to_cluster(cl, 'QxQStore', matrix(0, nrow=Q, ncol=Q))

}




