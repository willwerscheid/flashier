calc.obj <- function(flash, factor = NULL) {
  n.nonmissing <- get.n.nonmissing(flash)
  k            <- get.k(factor)

  KL <- sum(unlist(get.KL(flash)))
  if (!is.null(factor)) {
    KL <- KL + sum(unlist(get.KL(factor)))
    if (!is.new(factor))
      KL <- KL - sum(get.KLk(flash, k))
    est.tau <- get.est.tau(factor)
  } else {
    est.tau <- get.est.tau(flash)
  }

  # TODO: fixed S
  return(KL - 0.5 * sum(n.nonmissing * (log(2 * pi / est.tau) + 1)))
}

normal.means.loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) && s > 0
  x   <- x[idx]
  s   <- s[idx]
  Et  <- Et[idx]
  Et2 <- Et2[idx]
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}
