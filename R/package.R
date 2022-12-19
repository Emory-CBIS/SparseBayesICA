#' SparseBayesICA
#'
#' Package implementing codes for independent component analysis of multi-subject
#' fMRI data with sparse estimation of covariate effects.
#'
#' @docType package
#' @author Joshua Lukemire <joshua.lukemire@emory.edu>
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor kmeans prcomp rbeta rbinom rgamma rnorm runif sd var
#' @importFrom utils tail write.csv
#' @importFrom expm sqrtm
#' @importFrom reshape2 melt
#' @importFrom fastICA fastICA
#' @importFrom invgamma rinvgamma
#' @useDynLib SparseBayesICA
#' @name SparseBayesICA
#' @exportPattern "^[[:alpha:]]+"
NULL
