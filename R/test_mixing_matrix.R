if (0 == 1){


  #
  # Version 1: Hyperparameter Check ----
  #


  library(hcicaR)

  # data path
  dp <- "/Users/joshlukemire/Documents/Research/cbis/SparseBayesICA/simulation_demo_data"

  set.seed(100)




  N <- 10
  V <- 100000
  Q <- 3
  Htrue <- 5
  mu_h_true       <- c(2.0, 0.0, -3.0, 0.5, -0.5)
  sigma_h_sq_true <- c(2.0, 0.1, 0.5, 0.25, 0.8)
  CMtrue <- matrix(as.integer(sample(1:Htrue, size = Q*V, replace = TRUE, prob = c(100, 100, 30, 10, 5) )),
                   nrow = V,
                   ncol  = Q)

  # Generate true values for S0
  S0true <- matrix(0, nrow = V, ncol = Q)
  for (q in 1:Q){
    for (v in 1:V){
      S0true[v, q] <- rnorm(1, mean = mu_h_true[CMtrue[v, q]], sd = sqrt(sigma_h_sq_true[CMtrue[v, q]]) )
    }
  }

  P <- 2
  tau_sq_true    <- matrix(c(0.01, 0.001), nrow = Q, ncol = P)
  lambda_sq_true <- matrix(sample(c(0.01, 1000), size = Q * V * P, replace = TRUE, prob = c(0.9, 0.1)), nrow = V)
  # Generate the betas
  beta_true <- matrix(0, nrow = V, ncol = Q*P)
  cind <- 0
  for (q in 1:Q){
    for (p in 1:P){
      cind = cind + 1
      for (v in 1:V){
        beta_true[v,cind] <- rnorm(1, mean = 0, sd = sqrt(tau_sq_true[q,p] * lambda_sq_true[v, cind]) )
      }
    }
  }

  X <- matrix(rnorm(N*P), nrow = N)
  X_with_int  <- cbind(1, X)
  XtX         <- t(X_with_int) %*% X_with_int
  posterior_mean     <- rep(0, ncol(X_with_int))
  posterior_variance <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  prior_precision    <- matrix(0, ncol = ncol(X_with_int), nrow = ncol(X_with_int))
  draw <- rep(0, ncol(X_with_int))

  # Generate the "observed" data, this is directly A'Y
  AtY <- matrix(0, nrow = V, ncol = N * Q)
  Siall <- matrix(0, nrow = V, ncol = N * Q)
  inds <- 1:Q
  for (i in 1:N){
    Si <- S0true
    for (p in 1:P){
      bs <- seq(from = p, to = ncol(beta_true), by = P)
      Si <- Si + X[i, p] * beta_true[, bs]
    }
    Siall[,inds] <- Si
    AtY[,inds] <- Si + rnorm(Q*V)
    inds <- inds + Q
  }

  # Generate the mixing matrix for each subject
  Atrue <- NULL
  Ainit <- NULL
  for (i in 1:N){
    # random temporary mixing matrix
    Ai <- matrix(runif(Q^2, -3, 3), nrow = Q)

    # make orthonormal
    Ai <- Ai %*% solve(sqrtm(t(Ai) %*% Ai))

    # store
    Atrue <- rbind(Atrue, Ai)

    # random temporary mixing matrix
    Ai <- matrix(runif(Q^2, -3, 3), nrow = Q)

    # make orthonormal
    Ai <- Ai %*% solve(sqrtm(t(Ai) %*% Ai))

    # store
    Ainit <- rbind(Ainit, Ai)
  }

  inds <- 1:Q
  Y <- NULL
  for (i in 1:N){
    Yi <- AtY[, inds] %*% Atrue[inds, ]

    YYtemp <- t(Yi) %*% (Yi)
    ed <- eigen(YYtemp)
    Atrue[inds, ] <- t(ed$vectors) %*%  Atrue[inds, ]
    Yi <- Yi %*% ed$vectors

    Y <- cbind(Y, Yi)
    inds <- inds + Q
  }


  YYt <- matrix(0, nrow=N*Q, Q)
  YSiall      <- array(0, dim = c(Q*N, Q))
  for (iSubj in 1:N){
    inds <-  (Q*(iSubj-1)+1):(Q*iSubj)
    YYt[ inds, ] = t(Y[,inds]) %*% Y[,inds]
    YSiall[inds,] <- t(Y[,inds]) %*% (Siall[,inds ])
  }


  #
  # Sampler
  #

  Omega <- diag(rep(1, Q))
  YSiOmega<- array(0, dim = c(Q*N, Q))
  matmul_inplace(YSiOmega, YSiall, Omega)
  negYSiSigInvover2 <- -Omega/2
  # Angles for the mixing matrix update
  angle_range = seq(from = 0.0, to = 2*pi, length.out = 500)
  sin_theta <- sin(angle_range)
  cos_theta <- cos(angle_range)

  column_reduced_mixing_matrix <- matrix(0, nrow = Q, ncol = Q - 2)
  null_space <- matrix(0, nrow = Q, ncol = 2)

  gibbs_sample_mixing_matrix_inplace(Ainit,
                                     YSiOmega,
                                     negYSiSigInvover2,
                                     YYt,
                                     angle_range,
                                     sin_theta,
                                     cos_theta,
                                     column_reduced_mixing_matrix,
                                     null_space)

  Ainit[inds,] = rbmf.matrix.gibbs(A=YYt[inds,], B = negYSiSigInvover2, C = YSiOmega[inds,], Ainit[inds,])
  print( Ainit[inds,] )

cbind(Ainit, 0, Atrue)

# Least squares??
inds <- 1:Q
Als <- NULL
for (i in 1:N){
  proj <- Siall[,inds] %*% solve(t(Siall[,inds]) %*% Siall[,inds])
  au <- t(Y[,inds]) %*% proj
  aa <- au %*% sqrtm(solve(t(au)%*%au))
  Als <- rbind(Als, aa)
  inds <- inds + Q
}

cbind(Ainit, 0, Atrue, 0, Als)

inds <- 1:Q
t(Ainit[inds,]) %*% Atrue[inds,]
t(Als[inds,]) %*% Atrue[inds,]
inds <- inds + Q

}
