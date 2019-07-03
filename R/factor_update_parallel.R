update.factors.parallel <- function(flash, kset, cl) {
  for (n in 1:2) {
    ebnm.res <- solve.ebnm.parallel(n, flash, cl)

    flash <- set.EF(flash, sapply(ebnm.res, function(k) k$posterior$mean), n)
    flash <- set.EF2(flash, sapply(ebnm.res, function(k) k$posterior$second_moment), n)
    flash <- set.KL(flash, sapply(ebnm.res, function(k) k$KL), n)
    flash <- set.g(flash, lapply(ebnm.res, `[[`, "fitted_g"), n)
  }

  flash <- init.tau(flash) # need a faster way to do this; store Y2 and calc R^2 more efficiently
  flash <- set.obj(flash, calc.obj(flash))
  # TODO: set is.valid, is.zero

  # TODO: better convergence criteria?

  return(flash)
}

#' @importFrom parallel parLapply
#'
solve.ebnm.parallel <- function(n, flash, cl) {
  # TODO: add kset argument

  ebnm.args <- calc.all.ebnm.args(n, flash)

  ebnm.res <- parLapply(cl, ebnm.args, function(k) {
    # TODO: warmstarts?
    res    <- k$ebnm.fn(x = k$x, s = k$s, g = NULL, fixg = FALSE, output = k$output)
    res$KL <- res$log_likelihood - normal.means.loglik(k$x, k$s, res$posterior$mean,
                                                       res$posterior$second_moment)
    return(res)
  })

  return(ebnm.res)
}

calc.all.ebnm.args <- function(n, flash) {
  tau <- get.tau.for.ebnm.calc(flash)
  s2  <- calc.all.s2(n, flash, tau)
  x   <- calc.all.x(n, flash, s2, tau)
  s   <- sqrt(s2)

  # TODO: handle fixed factors, include.fixed option, handle exclusions

  return(lapply(1:ncol(x), function(k) {
    list(x = x[, k],
         s = s[, k],
         ebnm.fn = get.ebnm.fn.k(flash, k)[[n]],
         g = get.g.k(flash, k, n),
         output = default.output())
  }))
}

calc.all.s2 <- function(n, flash, tau) {
  other.n <- (1:2)[-n]
  if (is.tau.lowrank(flash)) {
    tau.EF2 <- as.vector(tau[[other.n]]) * get.EF2(flash, other.n)
    s2 <- 1 / outer(as.vector(tau[[n]]), colSums(tau.EF2))
  } else {
    # If tau is full-rank, then it has already been multiplied by Z:
    s2 <- 1 / nmode.prod.mat(tau, get.EF2(flash, other.n), other.n)
  }

  return(pmax(s2, 0))
}

calc.all.x <- function(n, flash, s2, tau) {
  Y <- get.Y(flash) # TODO: Implement for R (?)
  other.n <- (1:2)[-n]
  if (is.tau.lowrank(flash)) {
    EF.n       <- get.EF(flash, n)
    EF.other.n <- get.EF(flash, other.n)

    tau.EF <- as.vector(tau[[other.n]]) * EF.other.n
    Y.tau.EF <- nmode.prod.mat(Y, tau.EF, other.n)

    if (!any.missing(flash)) {
      EF.tau.EF  <- EF.n %*% (t(EF.other.n) %*% tau.EF)
      EFk.tau.EF <- t(t(EF.n) * colSums(EF.other.n * tau.EF))
    } else {
      # TODO: calc for missing data
      # ugly.mat <- matrix(apply(EF, 2, function(v) v * tau.EF), nrow = nrow(EF))
      # tmp <- nmode.prod.mat(Z, ugly.mat, other.n)
    }

    x <- as.vector(tau[[n]]) * (Y.tau.EF - EF.tau.EF + EFk.tau.EF)
  } else {
    # TODO: tau not lowrank
  }

  x <- s2 * x

  if (any(is.infinite(s2))) {
    x[is.infinite(s2)] <- 0
  }

  return(x)
}
