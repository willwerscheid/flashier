#' @importFrom parallel stopCluster
#'
update.factors.parallel <- function(flash, kset, cl) {
  is.zero <- is.zero(flash)
  is.valid <- is.valid(flash)

  for (n in 1:get.dim(flash)) {
    ebnm.res <- try(solve.ebnm.parallel(n, flash, kset, cl))
    if (inherits(ebnm.res, "try-error")) {
      parallel::stopCluster(cl)
      stop("Error encountered while backfitting mode ", n, " loadings.")
    }

    EF  <- get.EF(flash, n)
    EF2 <- get.EF2(flash, n)
    KL  <- get.KL(flash)[[n]]
    g   <- get.g(flash, n)

    EF[, kset]  <- sapply(ebnm.res, function(k) k$posterior$mean)
    EF2[, kset] <- sapply(ebnm.res, function(k) k$posterior$second_moment)
    KL[kset]    <- sapply(ebnm.res, function(k) k$KL)
    g[kset]     <- lapply(ebnm.res, function(k) k$fitted_g)

    flash <- set.EF(flash, EF, n)
    flash <- set.EF2(flash, EF2, n)
    flash <- set.KL(flash, KL, n)
    flash <- set.g(flash, g, n)

    is.zero[kset] <- (is.zero[kset]
                      | sapply(ebnm.res, function(k) all(k$posterior$mean == 0)))
  }

  is.valid[kset] <- TRUE

  flash <- set.is.zero(flash, is.zero)
  flash <- set.is.valid(flash, is.valid)

  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))

  return(flash)
}

#' @importFrom parallel parLapply
#'
solve.ebnm.parallel <- function(n, flash, kset, cl) {
  ebnm.args <- calc.all.ebnm.args(n, flash)
  ebnm.args <- ebnm.args[kset]
  ebnm.res  <- parallel::parLapply(cl, ebnm.args, parallel.ebnm.fn)

  return(ebnm.res)
}

calc.all.ebnm.args <- function(n, flash) {
  s2  <- calc.all.s2(n, flash)
  x   <- calc.all.x(n, flash, s2)
  s   <- sqrt(s2)

  return(lapply(1:ncol(x), function(k) {
    list(ebnm.fn = get.ebnm.fn.k(flash, k)[[n]],
         x = x[, k],
         s = s[k],
         g = get.g.k(flash, k, n),
         output = default.ebnm.output)
  }))
}

parallel.ebnm.fn <- function(k) {
  res    <- k$ebnm.fn(x = k$x, s = k$s, g_init = k$g, fix_g = FALSE, output = k$output)
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
