#' Fit All Available Diachronic Models
#'
#' Fit all currently available models to a dataset.
#'
#' @param data Dataset in long format
#' @param models Models to fit; by default, all available models are fit
#' @param reps Number of repetitions to take (with different initial parameter guesses)
#' @param sdfactor Initial parameter guess S.D. prefactor
#' @param method Optimization method to use
#' @param control Control parameters to optimization routine
#' @param mc.cores Number of parallel processor cores to use (Mac and Linux)
#' @return List with the following elements:
#' \describe{
#' \item{summary:}{Summary of fit; a dataframe combining the individual outputs for each model}
#' \item{coef:}{Maximum likelihood estimates of coefficients}
#' }
#' @details
#' For further explanation of the fitting routine, see `fit_dm`.
#'
#' @export
fit_dm_all <- function(data,
                       models = list("CRH", "VRH", "qCRH", "qVRH", "BRH"),
                       reps = 10,
                       sdfactor = 0.01,
                       method = "BFGS",
                       control = list(),
                       mc.cores = round(parallel::detectCores()/2)) {
  out <- parallel::mclapply(X=models, cre2::fit_dm, data=data, reps=reps, sdfactor=sdfactor, method=method, control=control, mc.cores=mc.cores)

  outs <- NULL
  outp <- vector("list", length(out))
  names(outp) <- models
  for (i in 1:length(out)) {
    outs <- rbind(outs, out[[i]]$summary)
    outp[[i]] <- out[[i]]$coef
  }

  list(summary=outs, coef=outp)
}


