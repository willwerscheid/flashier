#   Custom ebnm functions may also be used. They must accept parameters x (a
# vector of observations), s (a vector of standard errors), g (an
# initialization for the prior), fixg (a boolean indicating whether g should be
# considered fixed), and output. Normally, functions should return a named list
# containing elements postmean (the posterior mean for the observations),
# postmean2 (the posterior expected squared values), fitted_g (the fitted
# prior), and penloglik (the log likelihood for the observations given the
# prior). When output = "post_sampler", they should instead return a named list
# containing element post_sampler (a function that takes a single parameter
# nsamp and returns a matrix with nsamp rows and length(x) columns, with each
# row corresponding to a single sample from the posterior).

ebnm.ash = function(x, s, ash.param) {
  res <- do.call(ashr::ash,
                 c(list(betahat = as.vector(x),
                        sebetahat = as.vector(s)),
                   ash.param))

  if (!is.null(res$flash_data))
    res <- res$flash_data

  return(res)
}

ebnm.pn = function(x, s, ebnm.param) {
  res <- do.call(ebnm::ebnm,
                 c(list(x = as.vector(x),
                        s = as.vector(s)),
                   ebnm.param))

  if (!is.null(res$result)) {
    res$postmean  <- res$result$PosteriorMean
    res$postmean2 <- res$result$PosteriorMean2
    res$result    <- NULL
  }
  if (!is.null(res$loglik))
    names(res)[names(res) == "loglik"] <- "penloglik"

  return(res)
}
