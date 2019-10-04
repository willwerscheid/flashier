handle.prior.family <- function(prior.family, data.dim) {
  error.msg <- "Invalid argument to prior.family."

  if (!is.list(prior.family))
    stop(error.msg)

  # If a named list is found (rather than a list of named lists), it is
  #   interpreted as specifying the prior family (or families) for all factors.
  if (!is.null(names(prior.family[[1]])))
    prior.family <- list(prior.family)

  # Each top-level list element in prior.family corresponds to a factor.
  prior.family <- lapply(prior.family, function(k) {
    if (!is.list(k))
      stop(error.msg)
    if (length(k) == 1)
      k <- rep(k, data.dim)
    if (length(k) != data.dim)
      stop(error.msg)
    return(as.list(k))
  })

  prior.sign <- lapply(prior.family, lapply, `[[`, "sign")
  ebnm.fn    <- lapply(prior.family, lapply, `[[`, "ebnm.fn")

  return(list(prior.sign = prior.sign,
              ebnm.fn = ebnm.fn))
}
