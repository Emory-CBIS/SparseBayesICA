% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing_and_initialization.R
\name{obtain_s0_guess_fastica}
\alias{obtain_s0_guess_fastica}
\title{Use fastica to obtain a set of initial values for the population level components.}
\usage{
obtain_s0_guess_fastica(data)
}
\arguments{
\item{data}{A Q x V dataset that was obtained by performing PCA on the stacked,
dimension reduced subject level data. See \code{vignette("digits_example", package = "SparseBayesICA")}
for more details.}
}
\value{
A V x Q matrix, where each column contains a starting value
for an independent component.
}
\description{
`obtain_s0_guess_fastica()` takes a set of subject level time courses that
have been stacked across the time domain and reduced to dimension Q. It then
uses fastica to obtain the group level components, which can be used as an
initial guess for the population level components in the SparseBayes ICA model.
}
