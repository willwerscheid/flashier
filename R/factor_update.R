#' @exportS3Method NULL
update.factor <- function(factor, flash, update.tau = TRUE) {
  if (is.zero(factor))
    return(factor)

  for (n in 1:get.dim(flash))
    if (!is.zero(factor) && !all_fixed(factor, n))
      factor <- update.factor.one.n(factor, n, flash)

  if (update.tau)
    factor <- update.R2.tau.and.obj(factor, flash)

  return(factor)
}

#' @exportS3Method NULL
update.factor.one.n <- function(factor, n, flash) {
  nonmissing.thresh <- get.nonmissing.thresh(flash, n)
  if (nonmissing.thresh > 0) {
    prop.nonmissing <- calc.prop.nonmissing(factor, n, flash)
    # After exclusions, at least two loadings must remain.
    min.thresh <- sort(prop.nonmissing, decreasing = TRUE)[2]
    nonmissing.thresh <- min(nonmissing.thresh, min.thresh)
    exclude.idx <- which(prop.nonmissing < nonmissing.thresh)
    factor <- set.exclusions(factor, exclude.idx, n)
  }

  ebnm.res <- solve.ebnm(factor, n, flash)

  if (only.update.subset(factor, n, flash)) {
    idx.subset          <- get.idx.subset(factor)
    new.EF              <- get.EF(factor, n)
    new.EF[idx.subset]  <- ebnm.res$posterior$mean
    new.EF2             <- get.EF2(factor, n)
    new.EF2[idx.subset] <- ebnm.res$posterior$second_moment
  } else {
    new.EF  <- ebnm.res$posterior$mean
    new.EF2 <- ebnm.res$posterior$second_moment
  }

  factor <- set.EF(factor, new.EF, n)
  factor <- set.EF2(factor, new.EF2, n)
  factor <- set.KL(factor, ebnm.res$KL, n)
  factor <- set.g(factor, ebnm.res$fitted_g, n)

  if (all(abs(new.EF) < .Machine$double.eps)) {
    report.zero.factor(get.verbose.lvl(flash), get.k(factor))
    factor <- set.to.zero(factor)
  }

  return(factor)
}

#' @exportS3Method NULL
update.R2.tau.and.obj <- function(factor, flash) {
  factor <- update.tau(factor, flash)
  factor <- set.obj(factor, calc.obj(flash, factor))
  factor <- set.to.valid(factor)
  return(factor)
}

calc.prop.nonmissing <- function(factor, n, flash) {
  Z   <- get.nonmissing(flash)
  EF2 <- r1.square(r1.drop.dim(get.EF(factor), n))
  return(nmode.prod.r1(Z, EF2, n) / r1.sum(EF2))
}

#' @exportS3Method NULL
solve.ebnm <- function(factor, n, flash, output = default.ebnm.output) {
  fix.dim <- get.fix.dim(factor)
  if (use.subsetted.flash.data(factor, n))
    factor <- add.subset.data(factor, flash, fix.dim, get.idx.subset(factor))

  ebnm.fn    <- get.ebnm.fn(flash, factor, n)
  incl.fixed <- add.fixed.to.ebnm.args(factor, n, flash, output)
  ebnm.args  <- calc.ebnm.args(factor, n, flash, incl.fixed)

  prev.g <- get.g(factor, n)
  if (!identical(output, default.ebnm.output) && !is.null(prev.g)) {
    g    <- prev.g
    fixg <- TRUE
    ignored.warnings <- "mode and scale parameters are ignored"
  } else if (((is.new(factor) && warmstart.greedy(flash)) ||
              (!is.new(factor) && warmstart.backfits(flash)))
             && !is.null(prev.g)
             && warmstart.sanity.check(prev.g, ebnm.args$x, ebnm.args$s)) {
    g    <- prev.g
    fixg <- FALSE
    ignored.warnings <- NULL
  } else {
    g    <- NULL
    fixg <- FALSE
    ignored.warnings <- NULL
  }

  # This code block can be removed when ebnm is able to handle SEs equal to
  #   zero.
  if (identical(output, "lfsr")) {
    s.orig <- ebnm.args$s
    if (any(s.orig <= 0)) {
      ebnm.args$x <- ebnm.args$x[s.orig > 0]
      ebnm.args$s <- ebnm.args$s[s.orig > 0]
    }
  }

  ebnm.fn.ignore.warn <- function(g) {
    withCallingHandlers(
      ebnm.fn(
        x = ebnm.args$x,
        s = ebnm.args$s,
        g_init = g,
        fix_g = fixg,
        output = output
      ),
      warning = function(w) {
        if (!is.null(ignored.warnings)
            && any(startsWith(conditionMessage(w), ignored.warnings)))
          invokeRestart("muffleWarning")
      }
    )
  }

  # First attempt, possibly with warmstart.
  ebnm.res <- tryCatch(
    ebnm.fn.ignore.warn(g = g),
    error = function (cnd) {
      if (fixg | is.null(g)) {
        stop(paste("EBNM solver (no warmstart) threw error:", cnd))
      } else {
        NULL
      }
    }
  )

  # If warmstart failed, try again.
  if (is.null(ebnm.res)) {
    ebnm.res <- tryCatch(
      ebnm.fn.ignore.warn(g = NULL),
      error = function (cnd) {
        stop(paste("EBNM solver (after failed warmstart) threw error:", cnd))
      }
    )
  }

  if (identical(output, default.ebnm.output)) {
    ebnm.res$KL <- (ebnm.res$log_likelihood
                    - normal.means.loglik(ebnm.args$x,
                                          ebnm.args$s,
                                          ebnm.res$posterior$mean,
                                          ebnm.res$posterior$second_moment))
  }

  # This code block can also be removed after ebnm has been updated.
  if (identical(output, "lfsr")) {
    if (any(s.orig <= 0)) {
      lfsr <- rep(NA, length(s.orig))
      lfsr[s.orig > 0] <- ebnm.res$posterior$lfsr
      ebnm.res <- list(posterior = data.frame(lfsr = lfsr))
    }
  }

  return(ebnm.res)
}

calc.ebnm.args <- function(factor, n, flash, include.fixed) {
  tau <- get.tau.for.ebnm.calc(flash, tau = get.tau(factor))
  if (n %in% get.fix.dim(factor))
    tau <- full.or.lowrank.subset(tau, n, get.idx.subset(factor))

  s2 <- calc.s2(factor, n, flash, tau)
  x  <- calc.x(factor, n, flash, s2, tau)

  if (include.fixed) {
    idx.subset         <- get.idx.subset(factor)
    all.x              <- get.EF(factor, n)
    all.x[idx.subset]  <- x
    all.s2             <- rep(0, length(all.x))
    all.s2[idx.subset] <- s2
  } else {
    all.x  <- x
    all.s2 <- s2
  }

  all.s2[get.exclusions(factor, n)] <- Inf

  return(list(x = all.x, s = sqrt(all.s2)))
}

calc.s2 <- function(factor, n, flash, tau) {
  if (use.subsetted.flash.data(factor, n)) {
    Z <- get.Z.subset(factor)
  } else {
    Z <- get.nonmissing(flash)
  }

  factor.EF2 <- get.EF2(factor)
  if (n %in% get.fix.dim(factor))
    factor.EF2 <- r1.subset(factor.EF2, n, get.idx.subset(factor))

  if (is.tau.lowrank(flash)) {
    s2 <- 1 / premult.nmode.prod.r1(Z, tau, factor.EF2[-n], n)
  } else {
    # If tau is full-rank, then it has already been multiplied by Z:
    s2 <- 1 / nmode.prod.r1(tau, factor.EF2[-n], n)
  }

  return(pmax(s2, 0))
}

calc.x <- function(factor, n, flash, s2, tau) {
  if (use.subsetted.flash.data(factor, n)) {
    R        <- get.R.subset(factor)
    Y        <- get.Y.subset(factor)
    Z        <- get.Z.subset(factor)
    flash.EF <- get.EF.subset(factor)
  } else {
    R        <- get.R(flash)
    Y        <- get.Y(flash)
    Z        <- get.nonmissing(flash)
    flash.EF <- get.EF(flash)
  }

  k <- get.k(factor)
  if (uses.R(flash) && !is.new(factor)) {
    flash.EFk <- get.EF.k(flash, k)
    if (use.subsetted.flash.data(factor, n))
      flash.EFk <- r1.subset(flash.EFk, n, get.idx.subset(factor))
  }

  factor.EF <- get.EF(factor)
  if (n %in% get.fix.dim(factor))
    factor.EF <- r1.subset(factor.EF, n, get.idx.subset(factor))

  if (uses.R(flash)) {
    x <- premult.nmode.prod.r1(R, tau, factor.EF[-n], n)
    if (!is.new(factor)) {
      if (is.tau.lowrank(flash)) {
        EFk.tau <- elemwise.prod.lowrank.r1(tau, flash.EFk)
        x <- x + premult.nmode.prod.r1(Z, EFk.tau, factor.EF[-n], n)
      } else {
        EFk.tau <- elemwise.prod.fullrank.r1(tau, flash.EFk)
        x <- x + nmode.prod.r1(EFk.tau, factor.EF[-n], n)
      }
    }
  } else {
    x <- premult.nmode.prod.r1(Y, tau, factor.EF[-n], n)
    if (is.tau.lowrank(flash)) {
      flash.EF.tau <- lowranks.prod(tau, flash.EF, broadcast = TRUE)
      if (!is.new(factor))
        flash.EF.tau <- lowrank.drop.k(flash.EF.tau, k)
      x <- x - premult.nmode.prod.r1(Z, flash.EF.tau, factor.EF[-n], n)
    } else {
      flash.EF.minus.k <- flash.EF
      if (!is.new(factor))
        flash.EF.minus.k <- lowrank.drop.k(flash.EF, k)
      flash.EF.tau <- elemwise.prod.fullrank.lowrank(tau, flash.EF.minus.k)
      x <- x - nmode.prod.r1(flash.EF.tau, factor.EF[-n], n)
    }
  }

  x <- s2 * x

  if (any(is.infinite(s2))) {
    x[is.infinite(s2)] <- 0
  }

  return(x)
}

only.update.subset <- function(factor, n, flash) {
  return((n %in% get.fix.dim(factor)) && !use.fixed.to.est.g(flash))
}

use.subsetted.flash.data <- function(factor, n) {
  return((n %in% get.fix.dim(factor)) && !all_fixed(factor, n))
}

add.fixed.to.ebnm.args <- function(factor, n, flash, output) {
  return((n %in% get.fix.dim(factor))
          && (use.fixed.to.est.g(flash)
              || !identical(output, default.ebnm.output)))
}
