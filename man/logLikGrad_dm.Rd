% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ll.R
\name{logLikGrad_dm}
\alias{logLikGrad_dm}
\title{Gradients of Log-Likelihoods}
\usage{
logLikGrad_dm(x, data, model, inverse = TRUE)
}
\arguments{
\item{x}{Parameter vector}

\item{data}{Data set (in long format)}

\item{model}{Model}

\item{inverse}{Whether to return the additive inverse of the function value}
}
\description{
Compute the gradient of the log-likelihood (or its additive inverse) for diachronic models
at a given parameter value.
}
\details{
The option \code{inverse = TRUE} is specified by default, since most optimization algorithms
find the minimum rather than the maximum. For the possible values of the \code{model} argument,
see \code{fit_dm}.
}
\seealso{
\code{fit_dm}
}
