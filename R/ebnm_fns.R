#' Flashier EBNM functions
#'
#' Functions used to solve the empirical Bayes normal means problem.
#'   \code{ebnm.pn} is a wrapper to \code{ebnm::ebnm}. \code{ebnm.ash} calls
#'   \code{ashr::ash}.
#'
#' Custom EBNM functions may also be used. They must accept parameters
#'   \code{x}, \code{s}, \code{g} (with default \code{g = NULL}),
#'   \code{fixg} (with default \code{fixg = FALSE}), and \code{output}. When
#'   \code{output = "flash.data"}, they should result a list containing
#'   elements \code{fitted_g} (the fitted prior), \code{loglik} (the
#'   log-likelihood for the observations given the prior), and \code{result},
#'   which must itself be a list containing elements \code{PosteriorMean}
#'   (the posterior means for the observations) and \code{PosteriorMean2}
#'   (the posterior expected squared values). When \code{output = "sampler"},
#'   they should return a list containing element \code{post_sampler} (a
#'   function that takes a single parameter \code{nsamp} and returns a matrix
#'   with \code{nsamp} rows and \code{length(x)} columns, with each row
#'   corresponding to a single sample from the posterior). Finally, when
#'   \code{output = "lfsr"}, they should return a list containing element
#'   \code{result}, which is itself a list containing element \code{lfsr},
#'   a vector of local false sign rates. For both \code{output = "sampler"}
#'   and \code{output = "lfsr"}, a \code{NULL} value may be returned without
#'   otherwise affecting the fitting process.
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
#' @param output Can be set to \code{"flash.data"} to return summary data used
#'   to update factors; \code{"sampler"} to return a function that will sample
#'   from the posterior; or \code{"lfsr"} to return a vector of local false
#'   sign rates.
#'
#' @rdname ebnm.fn
#'
#' @export
#'
ebnm.pn = function(x, s, g = NULL, fixg = FALSE, output = "flash.data", ...) {
  output <- switch(output,
                   flash.data = c("result", "fitted_g", "loglik"),
                   sampler = "post_sampler",
                   lfsr = "lfsr")

  res <- ebnm::ebnm(x = x, s = s, g = g, fixg = fixg, output = output, ...)

  return(res)
}

#' @rdname ebnm.fn
#'
#' @export
#'
ebnm.ash = function(x, s, g = NULL, fixg = FALSE, output = "flash.data", ...) {
  output <- switch(output,
                   flash.data = c("PosteriorMean", "PosteriorSD",
                                  "fitted_g", "loglik"),
                   sampler = "post_sampler",
                   lfsr = "lfsr")

  res <- ashr::ash(betahat = as.vector(x), sebetahat = as.vector(s),
                   g = g, fixg = fixg, outputlevel = output, ...)

  if (!is.null(res$result)) {
    res$result$PosteriorMean2 <- (res$result$PosteriorMean^2
                                  + res$result$PosteriorSD^2)
  }

  return(res)
}
