#   Custom ebnm functions may also be used. They must accept parameters x (a
# vector of observations), s (a vector of standard errors), g (an
# initialization for the prior), fixg (a boolean indicating whether g should be
# considered fixed), and output. Normally, functions should return a named list
# containing elements postmean (the posterior mean for the observations),
# postmean2 (the posterior expected squared values), fitted_g (the fitted
# prior), and loglik (the log likelihood for the observations given the
# prior). When output = "sampler", they should instead return a named list
# containing element post_sampler (a function that takes a single parameter
# nsamp and returns a matrix with nsamp rows and length(x) columns, with each
# row corresponding to a single sample from the posterior).

ebnm.pn = function(x, s, g = NULL, fixg = FALSE, output = "flash.data", ...) {
  output <- switch(output,
                   flash.data = c("result", "fitted_g", "loglik"),
                   sampler = "post_sampler",
                   lfsr = "lfsr")

  res <- ebnm::ebnm(x = x, s = s, g = g, fixg = fixg, output = output, ...)

  return(res)
}


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
