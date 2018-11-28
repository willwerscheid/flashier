update.tau <- function(factor, flash) {
  if (is.tau.simple(flash)) {
    factor <- update.simple.tau(factor, flash)
  } else if (is.var.type.zero(flash)) {
    factor <- update.zero.tau(factor, flash)
  } else if (is.var.type.kronecker(flash)) {
    factor <- update.kronecker.tau(factor, flash)
  } else if (is.var.type.noisy(flash)) {
    factor <- update.noisy.tau(factor, flash)
  }
  return(factor)
}

update.simple.tau <- function(factor, flash) {
  delta.R2 <- calc.delta.R2.for.simple.tau(factor, flash)
  est.tau  <- estimate.simple.tau(flash, delta.R2)
  tau      <- tau.from.given.and.est(flash, est.tau)
  factor   <- set.delta.R2(factor, delta.R2)
  factor   <- set.est.tau(factor, est.tau)
  factor   <- set.tau(factor, tau)

  return(factor)
}

update.zero.tau <- function(factor, flash) {
  if (uses.R(flash))
    factor <- set.R(factor, calc.residuals(flash, factor))
  factor <- set.tau(factor, get.given.tau(flash))

  return(factor)
}

update.kronecker.tau <- function(factor, flash) {
  factor <- set.tau(factor, estimate.kronecker.tau(flash, factor))

  return(factor)
}

update.noisy.tau <- function(factor, flash) {
  if (uses.R(flash))
    factor <- set.R(factor, calc.residuals(flash, factor))
  noisy.tau <- estimate.noisy.tau(flash, factor)
  factor <- set.tau(factor, noisy.tau$tau)
  factor <- set.sum.tau.R2(factor, noisy.tau$sum.tau.R2)
}

estimate.simple.tau <- function(flash, delta.R2 = 0) {
  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

tau.from.given.and.est <- function(flash, est.tau) {
  given.tau <- get.given.tau(flash)

  if (is.var.type.zero(flash)) {
    # All variance is pre-specified ("zero" variance type):
    tau <- given.tau
  } else if (is.null(given.tau)) {
    # All variance is estimated:
    tau <- est.tau
  } else {
    # Otherwise both types of variance are used. In the simple case here, the
    #   total variance is estimated and the pre-specified variance functions
    #   as a minimum:
    tau <- pmin(given.tau, est.tau)
  }

  return(tau)
}

estimate.kronecker.tau <- function(flash, factor = NULL) {
  # Tol and maxiter are hardcoded (for now, at least).
  tol <- 1e-3
  maxiter <- 100

  tau <- get.tau(factor)
  if (is.null(tau))
    tau <- get.tau(flash)
  tau <- as.r1.tau(tau, flash)
  tau.dim <- get.est.tau.dim(flash)
  kron.nonmissing <- get.kron.nonmissing(flash)

  Z   <- get.nonmissing(flash)
  EF  <- get.new.EF(flash, factor)
  EF2 <- get.new.EF2(flash, factor)

  R2 <- get.latest.Rsquared(flash, factor, EF)

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.tau <- tau
    for (i in 1:length(tau.dim)) {
      n <- tau.dim[i]
      tau.R2 <- (nmode.prod.r1(R2, tau[-n], n)
                 + premult.nmode.prod.r1(Z, EF2, tau[-n], n)
                 - premult.nmode.prod.r1(Z, lowrank.square(EF), tau[-n], n))
      tau[[i]] <- kron.nonmissing[[i]] / tau.R2
    }
    max.chg <- calc.max.abs.chg(tau, old.tau)
  }

  return(tau)
}

estimate.noisy.tau <- function(flash, factor = NULL) {
  Z   <- get.nonmissing(flash)
  EF  <- get.new.EF(flash, factor)
  EF2 <- get.new.EF2(flash, factor)

  R2 <- get.latest.Rsquared(flash, factor, EF, set.missing.to.zero = FALSE)
  R2 <- R2 + lowrank.expand(EF2) - lowrank.expand(lowrank.square(EF))
  S2 <- get.given.S2(flash)

  var.type <- get.est.tau.dim(flash)
  any.missing <- !identical(Z, 1)

  if (var.type == 0) {
    if (any.missing) {
      R2.slice <- R2[as.logical(Z)]
      S2.slice <- S2[as.logical(Z)]
    } else {
      R2.slice <- R2
      S2.slice <- S2
    }
    est.S2 <- optimize.noisy(R2.slice, S2.slice)
    final.S2 <- S2 + est.S2
  } else {
    if (any.missing) {
      R2.slices <- R2
      S2.slices <- S2
      is.na(R2.slices) <- is.na(S2.slices) <- !Z
      R2.slices <- lapply(apply(R2.slices, var.type, list.with.no.NAs), unlist)
      S2.slices <- lapply(apply(S2.slices, var.type, list.with.no.NAs), unlist)
    } else {
      R2.slices <- lapply(apply(R2, var.type, list), unlist)
      S2.slices <- lapply(apply(S2, var.type, list), unlist)
    }
    est.S2 <- mapply(optimize.noisy, R2.slices, S2.slices)
    each <- prod(c(1, get.dims(flash)[1:get.dim(flash) < var.type]))
    final.S2 <- S2 + rep(est.S2, each = each)
  }

  tau <- Z / final.S2
  sum.tau.R2 <- sum(tau * R2)
  return(list(tau = tau, sum.tau.R2 = sum.tau.R2))
}

list.with.no.NAs <- function(x) list(x[!is.na(x)])

optimize.noisy <- function(R2, S2) {
  opt.res <- optimize(function(x) {sum(log(x + S2)) + sum(R2 / (x + S2))},
                      interval = c(0, mean(R2)))
  return(opt.res$minimum)
}

get.tau.for.ebnm.calc <- function(flash, tau = NULL) {
  if (is.null(tau))
    tau <- get.tau(flash)
  if (is(tau, "r1"))
    return(r1.to.lowrank(tau, flash))
  if (!is.tau.simple(flash))
    return(tau)
  return(as.lowrank.tau(tau, flash))
}

as.lowrank.tau <- function(tau, flash) {
  tau.n <- get.est.tau.dim(flash)
  tau.lowrank <- lapply(as.list(get.dims(flash)),
                        function(dim) {matrix(1, nrow = dim, ncol = 1)})
  if (is.null(tau.n) || (length(tau.n) == 1)) {
    n <- get.R2.n(flash)
    tau.lowrank[[n]] <- tau * tau.lowrank[[n]]
  } else {
    tau.lowrank[tau.n] <- mapply(`*`, tau, tau.lowrank[tau.n])
  }

  class(tau.lowrank) <- "lowrank"
  return(tau.lowrank)
}

as.r1.tau <- function(tau, flash) {
  r1 <- as.list(rep(1, get.dim(flash)))
  r1[get.est.tau.dim(flash)] <- tau

  class(r1) <- "r1"
  return(r1)
}
