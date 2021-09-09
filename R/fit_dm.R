#' Fit a diachronic model
#'
#' Fit a constant or variable rate effect model to data using a maximum likelihood method.
#'
#' @param data Data set in long format, with columns `time`, `context` and `response`
#' @param model Model to fit: see Details for options
#' @param reps Number of repetitions to take (with different initial parameter guesses)
#' @param sdfactor Initial parameter guess S.D. prefactor
#' @param method Optimization method to use
#' @param control Control parameters to optimization routine
#' @param mc.cores Number of parallel processor cores to use (Mac and Linux only)
#' @return Data frame with following components:
#' \describe{
#' \item{model:}{Model fitted}
#' \item{N:}{Number of data points}
#' \item{C:}{Number of contexts}
#' \item{K:}{Number of model parameters}
#' \item{method:}{Numerical optimization method employed}
#' \item{reps:}{Number of fitting repetitions}
#' \item{runtime:}{Runtime (wall-clock) of fitting}
#' \item{convergence:}{Convergence information returned by `optim`: `0` indicates convergence}
#' \item{logL:}{Maximum of log-likelihood function found}
#' }
#' @details
#' Currently available models are `"CRH"` and `"VRH"` for the classical Constant and Variable Rate Hypotheses (Kroch 1989), `"qCRH"` and `"qVRH"` for their quadratic extensions (Kallel 2007), and `"BRH"` for the Biased Rate Hypothesis (Kauhanen and Walkden 2018). See Kauhanen (in prep) for details of these, including calculations for the likelihood functions and their gradients.
#'
#' The maximum likelihood estimate is found using `optim`. By default, the BFGS algorithm is assumed (this seems to perform the best across the board), though other methods may be specified through the `method` (see the `optim` documentation for the available options). The behaviour of the optimization algorithm can be further controlled through use of the `control` parameter.
#'
#' An initial parameter guess is found by fitting a logistic model of the form `response~time` without regard for contexts. To avoid convergence to local optima, the optimization is repeated from `reps` such initial guesses, drawing each parameter from a normal distribution whose mean is the parameter found by the contextless logistic regression and whose standard deviation equals `sdfactor` times that mean; coefficient estimates for the quadratic terms in the qCRH and qVRH models, as well as those for the bias parameters in the BRH model, however, start with zero estimates as this seems to best guarantee convergence. The repetitions are conducted in parallel using `mclapply` if a multicore processor is present and forking is possible (i.e. on Linux and Mac systems). Increasing the value of `reps` will lessen the probability of idiosyncratic fits but will increase computation time.
#' @references
#' Kallel, Amel (2007) The loss of negative concord in Standard English: internal factors. *Language Variation and Change*, 19, 27--49.
#'
#' Kauhanen, Henri (in prep) Grammar competition, speaker models and rates of change: a critical reappraisal of the Constant Rate Hypothesis.
#' 
#' Kauhanen, Henri & Walkden, George (2018) Deriving the Constant Rate Effect. *Natural Language & Linguistic Theory*, 36, 483--521.
#'
#' Kroch, Anthony S. (1989) Reflexes of grammar in patterns of language change. *Language Variation and Change*, 1, 199--244.
#' @seealso
#' `optim`
#' @export
fit_dm <- function(data,
                   model,
                   reps = 10,
                   sdfactor = 0.01,
                   method = "BFGS",
                   control = list(),
                   mc.cores = round(parallel::detectCores()/2)) {
  dl <- rep(list(data), reps)
  runtime <- system.time(out <- parallel::mclapply(dl, fit_underthehood, model, method, control, sdfactor, mc.cores=mc.cores))[3]
  maxlogL <- out[[1]]$summary$logL
  whichmax <- 1
  for (i in 1:length(out)) {
    thislogL <- out[[i]]$summary$logL
    if (thislogL > maxlogL) {
      maxlogL <- thislogL
      whichmax <- i
    }
  }
  out <- out[[whichmax]]

  outs <- out$summary
  outs$method <- method
  outs$reps <- reps
  outs$runtime <- runtime
  out$summary <- outs[, c("model", "N", "C", "K", "method", "reps", "runtime", "convergence", "logL")]

  out
}


fit_underthehood <- function(data,
                             model,
                             method = "BFGS",
                             control = list(),
                             sdfactor = 0.01) {
  mod <- stats::glm(response~time, data, family=stats::binomial)
  slope <- mod$coefficients[2]
  intercept <- mod$coefficients[1]
  slope <- stats::rnorm(n=1, mean=slope, sd=sdfactor*abs(slope))
  intercept <- stats::rnorm(n=1, mean=intercept, sd=sdfactor*abs(intercept))
  qslope <- 0.0

  C <- length(unique(data$context))
  N <- nrow(data)

  if (model == "CRH") {
    mod <- stats::optim(par=c(slope, rep(intercept, C)), fn=logLik_dm, gr=logLikGrad_dm, model=model, data=data, method=method, control=control)
    m <- C + 1
  } else if (model == "qCRH") {
    mod <- stats::optim(par=c(qslope, slope, rep(intercept, C)), fn=logLik_dm, gr=logLikGrad_dm, model=model, data=data, method=method, control=control)
    m <- C + 2
  } else if (model == "BRH") {
    mod <- stats::optim(par=c(slope, intercept, stats::rnorm(n=C, mean=0, sd=0.0)), fn=logLik_dm, gr=logLikGrad_dm, model=model, data=data, method=method, control=control)
    m <- C + 2
  } else if (model == "VRH") {
    mod <- stats::optim(par=c(rep(slope, C), rep(intercept, C)), fn=logLik_dm, gr=logLikGrad_dm, model=model, data=data, method=method, control=control)
    m <- 2*C
  } else if (model == "qVRH") {
    mod <- stats::optim(par=c(rep(qslope, C), rep(slope, C), rep(intercept, C)), fn=logLik_dm, gr=logLikGrad_dm, model=model, data=data, method=method, control=control)
    m <- 3*C
  } else {
    stop("Invalid 'model' specified.")
  }

  logL <- -mod$value
  conv <- mod$convergence

  pars <- data.frame(t(mod$par))
  if (model == "CRH") {
    names(pars) <- c("s", paste0("k", 1:C))
  } else if (model == "qCRH") {
    names(pars) <- c("q", "s", paste0("k", 1:C))
  } else if (model == "BRH") {
    names(pars) <- c("s", "k", paste0("b", 1:C))
  } else if (model == "VRH") {
    names(pars) <- c(paste0("s", 1:C), paste0("k", 1:C))
  } else if (model == "qVRH") {
    names(pars) <- c(paste0("q", 1:C), paste0("s", 1:C), paste0("k", 1:C))
  }

  list(summary=data.frame(model=model, N=N, C=C, K=m, convergence=conv, logL=logL), parameters=pars)
}
