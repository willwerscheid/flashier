#' Wrapper to \code{ashr::ash}
#'
#' \code{flashier} is designed to use \code{ebnm::ebnm} to solve the empirical
#'   Bayes normal means problem. This function wraps \code{ashr::ash} so that
#'   it can use the \code{ebnm} interface.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard errors (or a scalar when all standard errors
#'   are the same).
#'
#' @param g An initialization for the prior.
#'
#' @param fixg A boolean indicating whether \code{g} should be considered
#'   fixed.
#'
#' @param output Can be set to \code{default.output()} to return summary data
#'   used to update factors; \code{"post_sampler"} to return a function that
#'   will sample from the posterior; or \code{"lfsr"} to return a vector of
#'   local false sign rates.
#'
#' @export
#'
ebnm.ash = function(x, s, g, fixg, output, ...) {
  if (identical(output, default.output())) {
    output <- c("PosteriorMean", "PosteriorSD", "fitted_g", "loglik")
  }

  res <- ashr::ash(betahat = as.vector(x), sebetahat = as.vector(s),
                   g = g, fixg = fixg, outputlevel = output, ...)

  if (!is.null(res$result) && !is.null(res$result$PosteriorMean)) {
    res$result$PosteriorMean2 <- (res$result$PosteriorMean^2
                                  + res$result$PosteriorSD^2)
  }

  return(res)
}
