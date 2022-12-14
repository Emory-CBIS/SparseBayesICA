# sample_spatial_maps_inplace(Eigen::Map<Eigen::MatrixXd> & S0,
#                             Eigen::Map<Eigen::MatrixXd> & beta,
#                             Eigen::Map<Eigen::VectorXd> & posterior_mean,
#                             Eigen::Map<Eigen::MatrixXd> & posterior_variance,
#                             Eigen::Map<Eigen::VectorXd> & draw,
#                             const Eigen::MatrixXd & AtY,
#                             const Eigen::MatrixXd & X,
#                             const Eigen::MatrixXd & XtX,
#                             const Eigen::VectorXd & sigma_sq_q,
#                             const Eigen::MatrixXd & tau_sq,
#                             const Eigen::MatrixXd & lambda_sq,
#                             const Eigen::VectorXd & mu_h_div_sigma_sq_h,
#                             const Eigen::VectorXd & sigma_sq_h,
#                             const Eigen::MatrixXi & cluster_membership)
#
#

library(hcicaR)

V = 4
Q = 3
P = 2
N = 3
S0 = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow=4) * 1.0
Beta = matrix( 1:(V*Q*P), nrow = V )  * 1.0
posterior_mean = rep(0, P + 1)  * 1.0
posterior_variance = matrix(0, nrow  =P+1, ncol = P+1)  * 1.0
draw = rep(0, P + 1)  * 1.0

AtY <- matrix(1:12, nrow=4, ncol = 9) * 1.0

AtY <-  matrix(c(1:12,
               c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*3,
               c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)*3), nrow=4, ncol = 9) * 1.0
#X   <- cbind(1, c(1, 2, 3), c(5, 6, 7))  * 1.0
X   <- cbind(1, c(0,1,1), c(0, -1, 1))  * 1.0
XtX <- t(X) %*% X  * 1.0
sigma_sq_q <- c(1.0, 1.0, 1.0)  * 0.01
tau_sq     <- matrix(100000.5, nrow = Q, ncol = P)
lambda_sq  <- matrix(1000000.1, nrow = V, ncol = P*Q)  * 1.0
mu_h_div_sigma_sq_h <- c(0.0, 0.0, 0.0)  * 1.0
sigma_sq_h <- c(100.0, 100.0, 100.0)  * 1.0
cluster_membership <- matrix(0:2, nrow=V, ncol = Q)  * 1.0
prior_precision <- matrix(0, nrow = 3, ncol = 3)

sample_spatial_maps_inplace(S0, Beta, posterior_mean, posterior_variance, prior_precision, draw,
                            AtY, X, XtX, sigma_sq_q, tau_sq, lambda_sq,
                            mu_h_div_sigma_sq_h, sigma_sq_h, (cluster_membership))
S0

Beta




test_that("mult", {

  sample_spatial_maps_inplace(S0, Beta, posterior_mean, posterior_variance, draw,
                              AtY, X, XtX, sigma_sq_q, tau_sq, lambda_sq,
                              mu_h_div_sigma_sq_h, sigma_sq_h, (cluster_membership))
  S0

  expect_equal(2 * 2, 4)
})



library(hcicaR)

xxx = load("/Users/joshlukemire/Desktop/debugState.Rdata")

sample_spatial_maps_inplace(S0, Beta, posterior_mean, posterior_variance, draw,
                            AtY, X_with_int, XtX, sigma_q_sq,
                            tau_sq, lambda_sq,
                            mu_h_div_sigma_sq_h, sigma_sq_h, cluster_memberships)

S0[1:4, 1:4]


