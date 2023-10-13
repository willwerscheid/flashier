calc.obj <- function(flash, factor = NULL) {
  KL <- sum(unlist(get.KL(flash)))

  if (!is.null(factor)) {
    KL <- KL + sum(unlist(get.KL(factor)))
    if (!is.new(factor))
      KL <- KL - sum(get.KL.k(flash, get.k(factor)))
    est.tau    <- get.est.tau(factor)
    tau        <- get.tau(factor)
    sum.tau.R2 <- get.sum.tau.R2(factor)
  } else {
    est.tau    <- get.est.tau(flash)
    tau        <- get.tau(flash)
    sum.tau.R2 <- get.sum.tau.R2(flash)
  }

  if (is.tau.simple(flash)) {
    n.nonmissing <- get.n.nonmissing(flash)
    obj <- KL - 0.5 * sum(n.nonmissing * (log(2 * pi)
                                          - log(tau) + tau / est.tau))
  } else if (is.var.type.zero(flash)) {
    obj <- KL - 0.5 * (get.log.2pi.s2(flash)
                       + sum(calc.tau.R2(flash, factor, get.R2.n(flash))))
  } else if (is.var.type.kronecker(flash)) {
    n.nonmissing <- get.kron.nonmissing(flash)
    obj <- KL - 0.5 * (sum(n.nonmissing[[1]]) * (log(2 * pi) + 1)
                       - sum(unlist(n.nonmissing) * log(unlist(tau))))
  } else if (is.var.type.noisy(flash) || is.var.type.noisy.kron(flash)) {
    if (any_missing(flash))
      tau <- tau[tau > 0]
    obj <- KL - 0.5 * (sum(log(2 * pi / tau)) + sum.tau.R2)
  } else {
    # This error should never occur:
    stop("The requested variance structure has not yet been implemented.")
  }

  return(obj)
}

normal.means.loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x   <- x[idx]
  s   <- s[idx]
  Et  <- Et[idx]
  Et2 <- Et2[idx]

  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}
