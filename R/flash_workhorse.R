#' A horse that works for \code{\link{flashier}}
#'
#' Caveat emptor.
#'
#' @inheritParams flashier
#'
#' @param prior.sign A vector or list of vectors indicating the sign(s) of
#'   the priors (-1 for nonpositive, 1 for nonnegative, and 0 otherwise). Only
#'   used when initializing factors.
#'
#' @param ebnm.fn A list of lists giving the functions to be used to solve the
#'   Empirical Bayes normal means problem when updating each factor.
#'
#' @param ebnm.param A list of lists giving the parameters to be passed to the
#'   functions in \code{ebnm.fn}.
#'
#' @param backfit.kset Which factors to backfit. The - operator can be used to
#'   instead specify which factors not to backfit. For example,
#'   \code{backfit.kset = -(1:2)} will backfit all factors but the first two.
#'
#' @param backfit.order How to determine the order in which to backfit factors.
#'   \code{"montaigne"} goes after the factor that promises to yield the
#'   largest increase in the variational lower bound. Il faut courir au plus
#'   pressé.
#'
#' @param warmstart.backfits Whether to use the current prior parameters to
#'   initialize the solution to the empirical Bayes normal means problem.
#'
#' @param backfit.after A vector of factor indices. A backfit will be
#'   performed each time one of these factors have been added.
#'
#' @param backfit.every After the last factor in \code{backfit.after} is added,
#'   an additional backfit will be performed after every \code{backfit.every}th
#'   factor is added.
#'
#' @param final.backfit Whether to perform a final backfit.
#'
#' @param nullchk.after Similar to \code{backfit.after}, but for nullchecks.
#'
#' @param nullchk.every Similar to \code{backfit.every}. For nullchecks.
#'
#' @param final.nullchk Whether to perform a final nullcheck.
#'
#' @param restart.after.nullchk Whether to continue fitting if the final
#'   nullcheck fails.
#'
#' @param conv.crit.fn The function to use to determine whether convergence has
#'   occurred. Used for both new factors and backfits.
#'
#' @param verbose.fns A vector of functions. Used to calculate values to
#'   output after each factor update.
#'
#' @param verbose.colnames A vector of column names, one for each function in
#'   \code{verbose.fns}.
#'
#' @param verbose.colwidths A vector of column widths.
#'
#' @param output.lvl What to include in the returned flash object. 0 = raw fit
#'   only; 1 = trimmed fit, no sampler; 2 trimmed fit and sampler; 3 = raw fit
#'   and sampler; 4 = raw fit, sampler, and lfsr (currently in beta).
#'
#' @param EF.init A list of matrices, one for each dimension. Each matrix
#'   should have k columns, one for each factor. New factors are initialized
#'   at these values.
#'
#' @param fix.dim A list of integers, one for each fixed factor. Specifies the
#'   dimension along which the factor is (partially) fixed.
#'
#' @param fix.idx A list of vectors, one for each fixed factor. Gives the
#'   indices of the loadings to fix.
#'
#' @param fix.vals A list of vectors that gives the values at which to fix
#'   the loadings specified by \code{fix.idx}.
#'
#' @param use.fixed.to.est.g Whether to include fixed values when estimating
#'   priors on (partially) fixed factors.
#'
#' @param nullchk.fixed.factors Whether to nullcheck fixed factors.
#'
#' @param init.fn The function to use to initialize factors.
#'
#' @param init.maxiter The maximum number of iterations when initializing
#'   factors.
#'
#' @param init.tol The maximum absolute change in normalized loadings that can
#'   occur before initialization is considered complete.
#'
#' @param greedy.maxiter The maximum number of iterations when optimizing a
#'   greedily added factor.
#'
#' @param greedy.tol The convergence tolerance when optimizing a greedily added
#'   factor.
#'
#' @param fixed.maxiter The maximum number of iterations when optimizing a
#'   newly added fixed factor.
#'
#' @param fixed.reltol The convergence tolerance (relative to greedy.tol) when
#'   optimizing a newly added fixed factor.
#'
#' @param backfit.maxiter Maximum iterations for final backfits.
#'
#' @param backfit.reltol Convergence tolerance for final backfits, relative to
#'   greedy.tol.
#'
#' @param inner.backfit.maxiter Maximum iterations for intermediary backfits.
#'
#' @param inner.backfit.reltol Convergence tolerance for intermediary backfits,
#'   relative to greedy.tol.
#'
#' @param nonmissing.thresh A vector of thresholds, one for each mode. Each
#'   threshold sets the (weighted) proportion of data that must be
#'   nonmissing in a given matrix or array slice in order to estimate the
#'   corresponding factor loading.
#'
#' @param seed By default, a seed is set for reproducibility. Set to
#'   \code{NULL} for a random seed.
#'
#' @param use.R Whether to maintain a matrix of residuals throughout the
#'   fitting process. This usually requires much more memory and seldom offers
#'   much improvement in runtime.
#'
flash.workhorse <- function(data = NULL,
                            flash.init = NULL,
                            var.type = 0,
                            prior.sign = NULL,
                            ebnm.fn = ebnm.pn,
                            ebnm.param = list(),
                            greedy.Kmax = 100,
                            backfit.kset = NULL,
                            backfit.order = c("dropout",
                                              "sequential",
                                              "random",
                                              "montaigne",
                                              "parallel"),
                            warmstart.backfits = TRUE,
                            backfit.after = NULL,
                            backfit.every = NULL,
                            final.backfit = FALSE,
                            nullchk.after = NULL,
                            nullchk.every = NULL,
                            final.nullchk = TRUE,
                            restart.after.nullchk = TRUE,
                            conv.crit.fn = calc.obj.diff,
                            verbose.lvl = 1,
                            verbose.fns = NULL,
                            verbose.colnames = NULL,
                            verbose.colwidths = NULL,
                            output.lvl = 3,
                            EF.init = NULL,
                            fix.dim = NULL,
                            fix.idx = NULL,
                            fix.vals = NULL,
                            use.fixed.to.est.g = FALSE,
                            nullchk.fixed.factors = FALSE,
                            init.fn = NULL,
                            init.maxiter = 100,
                            init.tol = 1e-2,
                            greedy.maxiter = 500,
                            greedy.tol = NULL,
                            fixed.maxiter = greedy.maxiter,
                            fixed.reltol = 1,
                            backfit.maxiter = 100,
                            backfit.reltol = 1,
                            inner.backfit.maxiter = backfit.maxiter,
                            inner.backfit.reltol = backfit.reltol,
                            nonmissing.thresh = NULL,
                            seed = 666,
                            use.R = FALSE) {
  set.seed(seed)
  backfit.order <- match.arg(backfit.order)

  ## data should be NULL when bypassing initialization.
  if (!is.null(data) && force.use.R(data, var.type)) {
    if (!missing(use.R) && !use.R)
      stop("R must be used with the requested variance structure.")
    use.R <- TRUE
  }
  if (!is.null(data) && inherits(data$Y, "lowrank") && use.R)
    stop("R cannot be used with low-rank Y.")

  announce.flash.init(verbose.lvl)
  if (!is.null(flash.init)) {
    if (missing(fix.dim))
      fix.dim <- flash.init$fix.dim
    if (missing(fix.idx))
      fix.idx <- flash.init$fix.idx
    if (missing(fix.vals))
      fix.vals <- flash.init$fix.vals
    if (missing(warmstart.backfits))
      warmstart.backfits <- flash.init$warmstart.backfits
    if (missing(use.fixed.to.est.g))
      use.fixed.to.est.g <- flash.init$use.fixed.to.est.g
    if (missing(nonmissing.thresh))
      nonmissing.thresh <- flash.init$nonmissing.thresh
  }
  flash <- init.flash(flash.init,
                      data = data,
                      EF.init = EF.init,
                      est.tau.dim = var.type,
                      dim.signs = prior.sign,
                      ebnm.fn = ebnm.fn,
                      ebnm.param = ebnm.param,
                      warmstart.backfits = warmstart.backfits,
                      fix.dim = fix.dim,
                      fix.idx = fix.idx,
                      fix.vals = fix.vals,
                      use.fixed.to.est.g = use.fixed.to.est.g,
                      nonmissing.thresh = nonmissing.thresh,
                      use.R = use.R)

  if (is.null(greedy.tol)) {
    greedy.tol <- set.default.tol(flash)
    report.tol.setting(verbose.lvl, greedy.tol)
  }

  total.factors.added <- 0
  max.factors.to.add  <- greedy.Kmax + get.n.fixed.to.add(flash)
  when.to.backfit <- as.Kset(backfit.after, backfit.every, max.factors.to.add)
  when.to.nullchk <- as.Kset(nullchk.after, nullchk.every, max.factors.to.add)

  if (verbose.lvl == -1)
    print.tab.delim.table.header(verbose.colnames)

  continue.looping <- TRUE
  is.converged     <- TRUE

  while (continue.looping) {
    continue.looping <- FALSE
    flash <- clear.flags(flash)

    continue.adding <- (total.factors.added < max.factors.to.add)
    if (continue.adding) {
      is.fixed <- is.next.fixed(flash)
      announce.add.factor(verbose.lvl, k = get.next.k(flash))

      announce.factor.init(verbose.lvl)
      factor <- init.factor(flash, init.fn, init.tol, init.maxiter)

      if (is.fixed) {
        maxiter <- fixed.maxiter
        tol     <- greedy.tol * fixed.reltol
      } else {
        maxiter <- greedy.maxiter
        tol     <- greedy.tol
      }

      if (maxiter > 0) {
        announce.factor.opt(verbose.lvl)
        print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                           backfit = FALSE)

        iter <- 0
        conv.crit <- Inf
        while (conv.crit > tol && iter < maxiter) {
          iter <- iter + 1

          old.f    <- factor
          factor   <- update.factor(factor, flash)
          obj.diff <- get.obj(factor) - get.obj(old.f)
          if (is.obj.valid(old.f) && obj.diff < 0
              && identical(get.exclusions(old.f), get.exclusions(factor))) {
            report.greedy.obj.decrease(verbose.lvl, obj.diff)
            factor <- old.f
            break
          }
          info <- calc.update.info(factor, old.f, conv.crit.fn, verbose.fns)
          conv.crit <- get.conv.crit(info)
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info,
                            get.next.k(flash), backfit = FALSE)
        }

        if (iter == maxiter)
          is.converged <- FALSE

        if (get.obj(factor) > get.obj(flash) + greedy.tol
            || !is.obj.valid(flash, factor)
            || is.fixed) {
          flash <- add.new.factor.to.flash(factor, flash)
        } else {
          flash <- set.greedy.fail.flag(flash)
        }
      } else if (is.fixed) {
        # Add fixed factors even when maxiter = 0.
        flash <- add.new.factor.to.flash(factor, flash)
      }

      if (greedy.failed(flash)) {
        continue.adding <- FALSE
      } else {
        continue.looping    <- TRUE
        total.factors.added <- total.factors.added + 1
      }

      report.add.factor.result(verbose.lvl, greedy.failed(flash),
                               get.obj(flash))
    }

    if (continue.adding) {
      do.backfit <- total.factors.added %in% when.to.backfit
      do.nullchk <- total.factors.added %in% when.to.nullchk
      maxiter    <- inner.backfit.maxiter
      tol        <- greedy.tol * inner.backfit.reltol
    } else {
      do.backfit <- final.backfit
      do.nullchk <- final.nullchk
      maxiter    <- backfit.maxiter
      tol        <- greedy.tol * backfit.reltol
    }

    if (do.backfit) {
      kset <- 1:get.n.factors(flash)
      conv.crit <- rep(Inf, get.n.factors(flash))
      if (!is.null(backfit.kset)) {
        # Remove any k that haven't been added yet.
        ksubset <- intersect(backfit.kset, c(kset, -kset))
        kset <- kset[ksubset]
        conv.crit[setdiff(1:get.n.factors(flash), kset)] <- 0
      }

      announce.backfit(verbose.lvl, n.factors = length(kset))
      print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                         backfit = TRUE)

      iter <- 0
      old.obj <- get.obj(flash)
      while (iter < maxiter && max(conv.crit) > tol) {
        is.converged <- TRUE
        iter <- iter + 1

        if (backfit.order == "random") {
          kset <- sample(kset)
        } else if (backfit.order == "montaigne") {
          # Il faut courir au plus pressé.
          last.max <- kset
          kset <- which.max(conv.crit)
          # Don't go after the same factor twice in a row.
          if (identical(last.max, kset)) {
            tmp <- conv.crit
            # Settle for second best.
            tmp[kset] <- 0
            if (max(tmp) > tol)
              kset <- which.max(tmp)
          }
        } else if (backfit.order == "dropout") {
          kset <- kset[conv.crit[kset] > tol]
        }

        if (backfit.order == "parallel") {
          old.f <- flash
          flash <- update.factors.parallel(flash, kset)
          info  <- calc.update.info(flash, old.f,
                                    conv.crit.fn, verbose.fns)
          conv.crit <- abs(get.conv.crit(info))
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info,
                            k = "all", backfit = TRUE)
        } else {
          for (k in kset) {
            old.f <- flash
            flash <- update.existing.factor(flash, k, iter, verbose.lvl)
            info  <- calc.update.info(flash, old.f,
                                      conv.crit.fn, verbose.fns, k)
            conv.crit[k] <- get.conv.crit(info)
            print.table.entry(verbose.lvl, verbose.colwidths, iter, info,
                              k = k, backfit = TRUE)
          }
        }
      }

      if (iter == maxiter)
        is.converged <- FALSE

      if (get.obj(flash) > old.obj) {
        report.backfit.complete(verbose.lvl, get.obj(flash))
      }
    }

    if (do.nullchk && get.n.factors(flash) > 0) {
      nullchk.kset <- 1:get.n.factors(flash)
      if (!nullchk.fixed.factors)
        nullchk.kset <- setdiff(nullchk.kset, which.k.fixed(flash))

      announce.nullchk(verbose.lvl, n.factors = length(nullchk.kset))

      for (k in nullchk.kset)
        flash <- nullcheck.factor(flash, k, verbose.lvl, greedy.tol)
      if (nullchk.failed(flash)) {
        if (restart.after.nullchk)
          continue.looping <- TRUE
      } else if (length(nullchk.kset) > 0) {
        report.nullchk.success(verbose.lvl)
      }
    }
  }

  if (!is.converged)
    warning("Flash fit has not converged. Try backfitting the returned fit, ",
            "setting backfit.tol and backfit.maxiter as needed.")

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl)

  report.completion(verbose.lvl)
  return(flash)
}

force.use.R <- function(data, var.type) {
  tmp <- list()

  tmp$est.tau.dim <- var.type
  tmp$given.tau.dim <- get.given.tau.dim(data)
  tmp$given.tau <- get.given.tau(data)
  tmp$given.S2 <- get.given.S2(data)

  return(is.var.type.zero(tmp) && !is.tau.simple(tmp))
}

as.Kset <- function(after, every, maxiter) {
  if (is.null(every))
    return(after)
  max.after <- max(c(0, after))
  n.every <- floor((maxiter - max.after) / every)
  every.iters <- max.after + every * (1:n.every)
  return(c(after, every.iters))
}
