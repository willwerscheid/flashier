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
  # ashr takes parameter "outputlevel", not "output".
  ash.param$outputlevel <- ash.param$output
  ash.param$output <- NULL

  res <- do.call(ashr::ash,
                 c(list(betahat = as.vector(x),
                        sebetahat = as.vector(s)),
                   ash.param))

  if (!is.null(res$flash_data))
    res <- res$flash_data
  else if (!is.null(res$result))
    res <- res$result

  return(res)
}

ebnm.pn = function(x, s, ebnm.param) {
  # TODO: update ebnm to return lfsr; update description above.
  if ("lfsr" %in% ebnm.param$output) {
    res <- list()
    res$lfsr <- rep(NA, length(x))
    return(res)
  }

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

# Check that an initial value of g has a fighting chance of fitting the data.
warmstart.sanity.check = function(g, x, s) {
  # Find rough limits for the region where g can have significantly positive
  #   density (pi can be ignored because it will be re-estimated).
  if (inherits(g,"unimix")) {
    upper.grange = max(g$b)
    lower.grange = min(g$a)
  } else if (inherits(g, "normalmix")) {
    # In the normal and halfnormal cases, use an anti-conservative range. It's
    #   better to re-estimate the grid than to use a bad one.
    upper.grange = 2 * max(g$sd)
    lower.grange = -upper.grange
  } else if (inherits(g, "tnormalmix")) {
    upper.grange = 2 * max(g$sd[is.infinite(g$b)])
    lower.grange = -2 * max(g$sd[is.infinite(g$a)])
  } else {
    # Skip the sanity check if g is not from a recognized ashr class.
    return(TRUE)
  }

  # Find the most outlying data points and calculate p-values conditional
  #   on the true values lying at the limits of g's range.
  worst.lower = min((x - lower.grange) / s)
  lower.p = pnorm(worst.lower)
  worst.upper = max((x - upper.grange) / s)
  upper.p = 1 - pnorm(worst.upper)

  # In the +uniform case, the worst lower point is ignored because we can't
  #   do anything about it.
  if (length(g$a) > 1 && min(g$a) == max(g$a)) {
    lower.p = 1
  }
  # And similarly for the -uniform case.
  if (length(g$b) > 1 && min(g$b) == max(g$b)) {
    upper.p = 1
  }

  worst.p = min(lower.p, upper.p)

  # The case where g is a single (null) component and the data comes from
  #   the null model should pass with 95% probability.
  n = length(x)
  passes.sanity.check = (worst.p > 1 - 0.95^(1/n))

  return(passes.sanity.check)
}
