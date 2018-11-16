calc.obj <- function(flash, factor = NULL) {
  n.nonmissing <- get.n.nonmissing(flash)
  k            <- get.k(factor)

  KL <- sum(unlist(get.KL(flash)))
  if (!is.null(factor)) {
    KL <- KL + sum(unlist(get.KL(factor)))
    # During backfitting, the old KL values need to be subtracted:
    if (!is.null(get.k(factor)))
      KL <- KL - sum(get.KLk(flash, k))
    est.tau <- get.est.tau(factor)
  } else {
    est.tau <- get.est.tau(flash)
  }

  return(KL - 0.5 * sum(n.nonmissing * (log(2 * pi / est.tau) + 1)))
  # TODO: fixed S
}

normal.means.loglik <- function(x, s, Et, Et2) {
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}
