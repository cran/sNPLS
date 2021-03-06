% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sNPLS_fit.R
\name{repeat_cv}
\alias{repeat_cv}
\title{Repeated cross-validation for sNPLS models}
\usage{
repeat_cv(
  X_npls,
  Y_npls,
  ncomp = 1:3,
  samples = 20,
  keepJ = NULL,
  keepK = NULL,
  nfold = 10,
  times = 30,
  parallel = TRUE,
  method = "sNPLS",
  ...
)
}
\arguments{
\item{X_npls}{A three-way array containing the predictors.}

\item{Y_npls}{A matrix containing the response.}

\item{ncomp}{A vector with the different number of components to test}

\item{samples}{Number of samples for performing random search in continuous thresholding}

\item{keepJ}{A vector with the different number of selected variables to test in discrete thresholding}

\item{keepK}{A vector with the different number of selected 'times' to test in discrete thresholding}

\item{nfold}{Number of folds for the cross-validation}

\item{times}{Number of repetitions of the cross-validation}

\item{parallel}{Should the computations be performed in parallel? Set up strategy first with \code{future::plan()}}

\item{method}{Select between sNPLS, sNPLS-SR or sNPLS-VIP}

\item{...}{Further arguments passed to cv_snpls}
}
\value{
A density plot with the results of the cross-validation and an (invisible) \code{data.frame} with these results
}
\description{
Performs repeated cross-validatiodn and represents results in a plot
}
