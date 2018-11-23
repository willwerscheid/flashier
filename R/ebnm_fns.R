ebnm.ash = function(x, s, ash.param) {
  res <- do.call(ashr::ash,
                 c(list(betahat = as.vector(x),
                        sebetahat = as.vector(s)),
                   ash.param))

  return(res[[1]])
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
