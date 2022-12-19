#' Main function for the SparseBayes ICA analysis
#'
#' @details
#' `fit_SparseBayesICA` carries out the posterior computation for the
#' SparseBayes ICA model.
#'
#' @param whitened_data A list containing the prewhitened time courses for each
#' subject. This can be generated using \code{\link{preprocess_subject}}
#' @param X A matrix containing the covariates to be included in the model. These
#' should already be in your preferred coding scheme
#' @param initial_values initial values for the SparseBayes ICA model parameters.
#' Can be obtained using \code{\link{obtain_initial_guess_sparsebayes}}
#' @param nburnin number of burnin iterations for the MCMC
#' @param nmcmc number of mcmc iterations after burnin
#' @param concentration The concentration parameter for the DPM. Larger numbers encourage
#' more clusters.
#' @param print_freq How often to print progress. Default is every 100 iterations
#' @param parclust NOT YET IMPLEMENTED IN PACKAGE. Cluster object for running
#' SparseBayes ICA in parallel.
#' @return A list containing the prewhitnened data and the eigenvalues of the
#' \itemize{
#'   \item Ai_posterior_mean - array containing the estimate for the mixing matrix
#'   for each subject
#'   \item S0_posterior_mean - matrix where each column is the posterior mean
#'   for one of the independent components
#'   \item Beta_posterior_mean - array containing the posterior mean of the
#'   covariate effect estimates at each voxel and component
#'   \item Beta_pseudo_pvalue -  array containing the pseudo-pvalue for each
#'   covariate effect at each voxel and component. This comes from examining
#'   how often the marginal credible interval overlaps with zero.
#'   \item sigma_sq_q_posterior_mean - posterior mean for the variance
#'   \item tau_sq_posterior_mean - posterior mean for the global shrinkage parameter
#'   of the horseshoe
#'   \item lambda_sq_posterior_mean - posterior mean for the local shrinkage parameters
#'   of the horseshoe
#'   \item timing - tibble with timing information
#' }
#' @export
fit_SparseBayesICA <- function(whitened_data, X,
                               initial_values = NULL,
                                 nburnin = 5000, nmcmc = 5000,
                                 concentration = 1.0,
                                 print_freq = 100,
                                 parclust = NULL){

  # Check if a random initialization needs to be constructed
  if (is.null(initial_values)){
    initial_values <- random_initialization(whitened_data, X)
  }

  V <- dim(initial_values$S0)[1]
  Q <- dim(initial_values$S0)[2]

  # Reformat Y if needed
  if (is.list(whitened_data)){
    Y = array(0, dim = c(V, N*Q))
    qSeq <- 1:Q
    for (i in 1:N){
      Y[,qSeq] = whitened_data[[i]]
      qSeq <- qSeq + Q
    }
  }

  # Check if running in parallel
  if (is.null(parclust)){
    result = SparseBayesICASerial(Y, X, initial_values,
                          nburnin = nburnin, nmcmc = nmcmc,
                          concentration = concentration,
                          print_freq = print_freq)

  } else {
    stop("Parallel SparseBayes is WIP")
  }

  return(result)

}
