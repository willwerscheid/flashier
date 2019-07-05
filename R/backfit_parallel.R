#' @importFrom parallel stopCluster
#'
update.factors.parallel <- function(flash, kset, cl) {
  # TODO: check that this works for tensors, sparse and lowrank
  for (n in 1:get.dim(flash)) {
    ebnm.res <- try(solve.ebnm.parallel(n, flash, kset, cl))
    if (inherits(ebnm.res, "try-error")) {
      parallel::stopCluster(cl)
      stop("Error encountered while backfitting mode ", n, " loadings.")
    }

    flash <- set.EF(flash, sapply(ebnm.res, function(k) k$posterior$mean), n)
    flash <- set.EF2(flash, sapply(ebnm.res, function(k) k$posterior$second_moment), n)
    flash <- set.KL(flash, sapply(ebnm.res, function(k) k$KL), n)
    flash <- set.g(flash, lapply(ebnm.res, `[[`, "fitted_g"), n)
  }

  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))
  # TODO: set is.valid, is.zero

  # TODO: better convergence criteria?

  return(flash)
}

#' @importFrom parallel parLapply
#'
solve.ebnm.parallel <- function(n, flash, kset, cl) {
  # TODO: add kset argument

  ebnm.args <- calc.all.ebnm.args(n, flash)
  ebnm.res  <- parallel::parLapply(cl, ebnm.args, parallel.ebnm.fn)

  return(ebnm.res)
}

calc.all.ebnm.args <- function(n, flash) {
  s2  <- calc.all.s2(n, flash)
  x   <- calc.all.x(n, flash, s2)
  s   <- sqrt(s2)

  # TODO: handle fixed factors, include.fixed option, handle exclusions

  return(lapply(1:ncol(x), function(k) {
    list(ebnm.fn = get.ebnm.fn.k(flash, k)[[n]],
         x = x[, k],
         s = s[k],
         g = get.g.k(flash, k, n),
         output = default.output())
  }))
}

parallel.ebnm.fn <- function(k) {
  res    <- k$ebnm.fn(x = k$x, s = k$s, g = k$g, fixg = FALSE, output = k$output)
  res$KL <- res$log_likelihood - normal.means.loglik(k$x,
                                                     k$s,
                                                     res$posterior$mean,
                                                     res$posterior$second_moment)
  return(res)
}

# Returns a k-vector.
calc.all.s2 <- function(n, flash) {
  tau <- get.tau(flash)
  EF2 <- get.EF2(flash)

  EF2.sums <- Reduce(`*`, lapply(EF2[-n], colSums))
  s2       <- 1 / (tau * EF2.sums)

  return(pmax(s2, 0))
}

calc.all.x <- function(n, flash, s2) {
  Y   <- get.Y(flash)
  tau <- get.tau(flash)
  EF  <- get.EF(flash)

  EF.crossprod <- Reduce(`*`, lapply(EF[-n], crossprod))
  EF.colsums2  <- Reduce(`*`, lapply(EF[-n], function(x) colSums(x^2)))

  x <- nmode.prod.lowrank(Y, EF[-n], n)
  x <- x - EF[[n]] %*% EF.crossprod
  x <- x + EF[[n]] * matrix(EF.colsums2,
                            nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  x <- tau * x * matrix(s2, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  if (any(is.infinite(s2))) {
    x[is.infinite(s2)] <- 0
  }

  return(x)
}
