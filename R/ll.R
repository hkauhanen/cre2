#' Log-Likelihoods
#'
#' Compute the log-likelihood (or its additive inverse) for diachronic models
#' at a given parameter value.
#'
#' @param x Parameter vector
#' @param data Data set (in long format)
#' @param model Model
#' @param inverse Whether to return the additive inverse of the function value
#' @details
#' The option `inverse = TRUE` is specified by default, since most optimization algorithms
#' find the minimum rather than the maximum. For the possible values of the `model` argument,
#' see `fit_dm`.
#' @seealso
#' `fit_dm`
#'
#' @export
logLik_dm <- function(x,
                      data,
                      model,
                      inverse = TRUE) {
  # turn contexts into integers
  data$context <- as.integer(factor(data$context))

  # total number of contexts
  C <- max(data$context)

  if (model == "CRH") {
    FUN <- ll_classical
  } else if (model == "BRH") {
    FUN <- ll_bias
  } else if (model == "VRH") {
    FUN <- ll_VRE
  } else if (model == "qCRH") {
    FUN <- ll_qCRE
  } else if (model == "qVRH") {
    FUN <- ll_qVRE
  } else {
    stop("Invalid 'model' specified")
  }

  out <- sum(apply(data, FUN=FUN, MARGIN=1, theta=x, C=C))

  if (inverse) {
    out <- -out
  }

  out
}


#' Gradients of Log-Likelihoods
#'
#' Compute the gradient of the log-likelihood (or its additive inverse) for diachronic models
#' at a given parameter value.
#'
#' @param x Parameter vector
#' @param data Data set (in long format)
#' @param model Model
#' @param inverse Whether to return the additive inverse of the function value
#' @details
#' The option `inverse = TRUE` is specified by default, since most optimization algorithms
#' find the minimum rather than the maximum. For the possible values of the `model` argument,
#' see `fit_dm`.
#' @seealso
#' `fit_dm`
#'
#' @export
logLikGrad_dm <- function(x,
                          data,
                          model,
                          inverse = TRUE) {
  # turn contexts into integers
  data$context <- as.integer(factor(data$context))

  # total number of contexts
  C <- max(data$context)

  if (model == "CRH") {
    FUN <- gr_classical
  } else if (model == "BRH") {
    FUN <- gr_bias
  } else if (model == "VRH") {
    FUN <- gr_VRE
  } else if (model == "qCRH") {
    FUN <- gr_qCRE
  } else if (model == "qVRH") {
    FUN <- gr_qVRE
  } else {
    stop("Invalid 'model' specified")
  }

  out <- rowSums(apply(data, FUN=FUN, MARGIN=1, theta=x, C=C))

  if (inverse) {
    out <- -out
  }

  out
}


ll_classical <- function(X,
                         C,
                         theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[1]*X$time + theta[1 + X$context])))
  X$response*log(p) + (1 - X$response)*log(1 - p)
}


ll_qCRE <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[1]*X$time^2 + theta[2]*X$time + theta[2 + X$context])))
  X$response*log(p) + (1 - X$response)*log(1 - p)
}


ll_VRE <- function(X,
                   C,
                   theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[X$context]*X$time + theta[C + X$context])))
  X$response*log(p) + (1 - X$response)*log(1 - p)
}


ll_qVRE <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[X$context]*X$time^2 + theta[C + X$context]*X$time + theta[2*C + X$context])))
  X$response*log(p) + (1 - X$response)*log(1 - p)
}


ll_bias <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  sigma <- -1 + 2/(1 + exp(-theta[2 + X$context]))
  lambda <- 1/(1 + exp(-(theta[1]*X$time + theta[2])))
  p <- lambda + sigma*lambda*(1-lambda)
  X$response*log(p) + (1 - X$response)*log(1 - p)
}


gr_classical <- function(X,
                         C,
                         theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[1]*X$time + theta[1 + X$context])))
  slope <- (X$response*(1 - p) - (1 - X$response)*p)*X$time
  int <- rep(0, C)
  int[X$context] <- X$response*(1 - p) - (1 - X$response)*p
  c(slope, int)
}


gr_qCRE <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[1]*X$time^2 + theta[2]*X$time + theta[2 + X$context])))
  qslope <- (X$response*(1 - p) - (1 - X$response)*p)*X$time^2
  slope <- (X$response*(1 - p) - (1 - X$response)*p)*X$time
  int <- rep(0, C)
  int[X$context] <- X$response*(1 - p) - (1 - X$response)*p
  c(qslope, slope, int)
}


gr_VRE <- function(X,
                   C,
                   theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[X$context]*X$time + theta[C + X$context])))
  slope <- rep(0, C)
  slope[X$context] <- (X$response*(1 - p) - (1 - X$response)*p)*X$time
  int <- rep(0, C)
  int[X$context] <- X$response*(1 - p) - (1 - X$response)*p
  c(slope, int)
}


gr_qVRE <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  p <- 1/(1 + exp(-(theta[X$context]*X$time^2 + theta[C + X$context]*X$time + theta[2*C + X$context])))
  slope <- rep(0, C)
  slope[X$context] <- (X$response*(1 - p) - (1 - X$response)*p)*X$time
  qslope <- rep(0, C)
  qslope[X$context] <- (X$response*(1 - p) - (1 - X$response)*p)*X$time^2
  int <- rep(0, C)
  int[X$context] <- X$response*(1 - p) - (1 - X$response)*p
  c(qslope, slope, int)
}


gr_bias <- function(X,
                    C,
                    theta) {
  X <- as.data.frame(t(X))
  sigma <- -1 + 2/(1 + exp(-theta[2 + X$context]))
  lambda <- 1/(1 + exp(-(theta[1]*X$time + theta[2])))
  p <- lambda + sigma*lambda*(1-lambda)
  S <- 1 + sigma*(1 - 2*lambda)
  slope <- (X$response/p - (1 - X$response)/(1 - p))*lambda*(1 - lambda)*S*X$time
  int <- (X$response/p - (1 - X$response)/(1 - p))*lambda*(1 - lambda)*S
  bias <- rep(0, C)
  bias[X$context] <- (X$response/p - (1 - X$response)/(1 - p))*lambda*(1 - lambda)*((2*exp(-theta[2 + X$context]))/(1 + exp(-theta[2 + X$context]))^2)
  c(slope, int, bias)
}

