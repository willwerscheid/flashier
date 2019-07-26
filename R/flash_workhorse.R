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
#' @param extrapolate.greedy Whether to use an extrapolation technique
#'   inspired by Ang and Gillis (2019) to accelerate the fitting of greedy
#'   and fixed factors.
#'
#' @param backfit.kset Which factors to backfit. The - operator can be used to
#'   instead specify which factors not to backfit. For example,
#'   \code{backfit.kset = -(1:2)} will backfit all factors but the first two.
#'
#' @param backfit.method \code{"sequential"} updates each factor in order, one
#'   at a time. \code{"extrapolate"} does the same but uses an acceleration
#'   technique inspired by Ang and Gillis (2019). The extrapolation parameters
#'   can be set via \code{extrapolate.control}. \code{"dropout"} is similar to
#'   \code{"sequential"}, but stops updating individual factors once they are
#'   no longer changing very much. \code{"montaigne"} goes after the factor that
#'   promises to yield the largest increase in the variational lower bound.
#'   Il faut courir au plus pressé. \code{"random"} updates each factor one at
#'   a time, but re-orders them randomly after each backfit iteration.
#'   \code{"parallel"} does a simultaneous update of all factors. Unlike other
#'   methods, parallel backfits are not guaranteed to yield monotonic increases
#'   in the variational lower bound.
#'   The number of cores used by \code{"parallel"} can be set via the
#'   command \code{options("cl.cores", n.cores)}. The type of multicore
#'   cluster can be set via \code{options("cl.type", type)}. Typically,
#'   \code{cl.type = "FORK"} is more efficient on Unix-likes.
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
#' @param tol The convergence tolerance.
#'
#' @param fixed.maxiter The maximum number of iterations when optimizing a
#'   newly added fixed factor.
#'
#' @param fixed.reltol The convergence tolerance (relative to tol) when
#'   optimizing a newly added fixed factor.
#'
#' @param backfit.maxiter Maximum iterations for final backfits.
#'
#' @param backfit.reltol Convergence tolerance for final backfits, relative to
#'   tol.
#'
#' @param inner.backfit.maxiter Maximum iterations for intermediary backfits.
#'
#' @param inner.backfit.reltol Convergence tolerance for intermediary backfits,
#'   relative to tol.
#'
#' @param nullchk.reltol Tolerance for nullchecks, relative to tol.
#'
#' @param extrapolate.control Extrapolation parameters. TODO: describe.
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
#' @importFrom parallel makeCluster stopCluster
#'
flash.workhorse <- function(data = NULL,
                            init = NULL,
                            var.type = 0,
                            prior.sign = NULL,
                            ebnm.fn = ebnm::ebnm,
                            greedy.Kmax = 100,
                            extrapolate.greedy = TRUE,
                            backfit.kset = NULL,
                            backfit.method = c("extrapolate",
                                               "sequential",
                                               "dropout",
                                               "montaigne",
                                               "random",
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
                            tol = NULL,
                            fixed.maxiter = greedy.maxiter,
                            fixed.reltol = 1,
                            backfit.maxiter = 100,
                            backfit.reltol = 1,
                            inner.backfit.maxiter = backfit.maxiter,
                            inner.backfit.reltol = backfit.reltol,
                            nullchk.reltol = 1,
                            extrapolate.control = list(),
                            nonmissing.thresh = NULL,
                            seed = 666,
                            use.R = FALSE) {
  set.seed(seed)
  backfit.method <- match.arg(backfit.method)

  # data should be NULL when bypassing initialization.
  if (!is.null(data) && force.use.R(data, var.type)) {
    if (!missing(use.R) && !use.R)
      stop("R must be used with the requested variance structure.")
    use.R <- TRUE
  }
  if (!is.null(data) && inherits(data$Y, "lowrank") && use.R)
    stop("R cannot be used with low-rank Y.")

  announce.flash.init(verbose.lvl)
  if (!is.null(init)) {
    if (missing(fix.dim))
      fix.dim <- init$fix.dim
    if (missing(fix.idx))
      fix.idx <- init$fix.idx
    if (missing(fix.vals))
      fix.vals <- init$fix.vals
    if (missing(warmstart.backfits))
      warmstart.backfits <- init$warmstart.backfits
    if (missing(use.fixed.to.est.g))
      use.fixed.to.est.g <- init$use.fixed.to.est.g
    if (missing(nonmissing.thresh))
      nonmissing.thresh <- init$nonmissing.thresh
  }
  flash <- init.flash(init,
                      data = data,
                      EF.init = EF.init,
                      est.tau.dim = var.type,
                      dim.signs = prior.sign,
                      ebnm.fn = ebnm.fn,
                      warmstart.backfits = warmstart.backfits,
                      fix.dim = fix.dim,
                      fix.idx = fix.idx,
                      fix.vals = fix.vals,
                      use.fixed.to.est.g = use.fixed.to.est.g,
                      nonmissing.thresh = nonmissing.thresh,
                      use.R = use.R)

  # For parallel backfits to be worthwhile, there must be an efficient way to
  #   calculate R2 (that is, more efficient than explicitly forming the n x p
  #   matrix of residuals). When the estimated variance is not constant across
  #   all entries or if there is any missing data, then (at minimum) a n x k^2
  #   or p x k^2 matrix must be created, so it must be true that
  #   k^2 << max(n, p) for parallelization to yield any real benefits. But for
  #   small k, serial backfits are fast and yield monotonic increases in the
  #   ELBO, so parallelization is best avoided.
  if (backfit.method == "parallel") {
    if (get.dim(flash) > 2) {
      stop("Parallel backfits have not yet been implemented for tensors.")
    }
    if (any.missing(flash)) {
      stop("Parallel backfits have not been implemented for missing data.")
    }
    if (!store.R2.as.scalar(flash)) {
      stop("Parallel backfits can only be performed when the estimated ",
           "variance is constant across all entries.")
    }
    if (use.R) {
      stop("Parallel backfits require use.R = FALSE.")
    }
    if (length(which.k.fixed(flash)) > 0) {
      if (length(intersect(which.k.fixed(flash), backfit.kset)) > 0) {
        stop("Parallel backfits have not yet been implemented for fixed ",
             "factors.")
      }
      if (is.null(backfit.kset)) {
        warning("Parallel backfits have not yet been implemented for fixed ",
                "factors, so they will be removed from backfit.kset.")
      }
    }
  }

  if (is.null(tol)) {
    tol <- set.default.tol(flash)
    report.tol.setting(verbose.lvl, tol)
  }

  total.factors.added <- 0
  max.factors.to.add  <- greedy.Kmax + get.n.fixed.to.add(flash)
  when.to.backfit <- as.Kset(backfit.after, backfit.every, max.factors.to.add)
  when.to.nullchk <- as.Kset(nullchk.after, nullchk.every, max.factors.to.add)

  extrapolate.param <- set.extrapolate.param(extrapolate.control)

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
        add.tol <- tol * fixed.reltol
      } else {
        maxiter <- greedy.maxiter
        add.tol <- tol
      }

      if (maxiter > 0) {
        announce.factor.opt(verbose.lvl)
        print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                           backfit = FALSE)

        iter <- 0
        conv.crit <- Inf
        if (extrapolate.greedy) {
          old.f <- factor
          extrapolate.param <- init.beta(extrapolate.param)
        }
        while (conv.crit > add.tol && iter < maxiter) {
          iter <- iter + 1

          if (extrapolate.greedy) {
            proposed.factor <- extrapolate.f(factor, old.f, extrapolate.param)
            proposed.factor <- update.factor(proposed.factor, flash)
          }

          old.f <- factor
          if (!extrapolate.greedy) {
            factor <- update.factor(factor, flash)
          } else if (get.obj(proposed.factor) - get.obj(factor) < tol) {
            factor <- update.factor(factor, flash)
            extrapolate.param <- decelerate(extrapolate.param)
          } else {
            factor <- proposed.factor
            extrapolate.param <- accelerate(extrapolate.param)
          }

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

        if (get.obj(factor) > get.obj(flash) + tol
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
      do.backfit  <- total.factors.added %in% when.to.backfit
      do.nullchk  <- total.factors.added %in% when.to.nullchk
      maxiter     <- inner.backfit.maxiter
      backfit.tol <- tol * inner.backfit.reltol
    } else {
      do.backfit  <- final.backfit
      do.nullchk  <- final.nullchk
      maxiter     <- backfit.maxiter
      backfit.tol <- tol * backfit.reltol
    }

    if (do.backfit) {
      kset <- 1:get.n.factors(flash)
      conv.crit <- rep(Inf, get.n.factors(flash))
      if (!is.null(backfit.kset)) {
        if (is.function(backfit.kset)) {
          next.kset <- backfit.kset(get.n.factors(flash))
          # Remove any k that haven't been added yet.
          ksubset <- intersect(next.kset, c(kset, -kset))
        } else {
          ksubset <- intersect(backfit.kset, c(kset, -kset))
        }
        kset <- kset[ksubset]
        conv.crit[setdiff(1:get.n.factors(flash), kset)] <- 0
      }

      announce.backfit(verbose.lvl, n.factors = length(kset), backfit.tol)
      print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                         backfit = TRUE)

      if (backfit.method == "parallel") {
        # Remove zero factors and fixed factors.
        kset <- setdiff(kset, which(is.zero(flash)))
        kset <- setdiff(kset, which.k.fixed(flash))

        cl <- parallel::makeCluster(getOption("cl.cores", 2L),
                                    type = getOption("cl.type", "PSOCK"),
                                    useXDR = FALSE)
      }

      iter <- 0
      old.obj <- get.obj(flash)
      next.tol.target <- NULL
      if (backfit.method == "extrapolate") {
        extrapolate.param <- init.beta(extrapolate.param)
        old.f <- flash
      }
      while (iter < maxiter && max(conv.crit) > backfit.tol) {
        is.converged <- TRUE
        iter <- iter + 1

        if (backfit.method == "random") {
          kset <- sample(kset)
        } else if (backfit.method == "montaigne") {
          # Il faut courir au plus pressé.
          last.max <- kset
          kset <- which.max(conv.crit)
          # Don't go after the same factor twice in a row.
          if (identical(last.max, kset)) {
            tmp <- conv.crit
            # Settle for second best.
            tmp[kset] <- 0
            if (max(tmp) > backfit.tol)
              kset <- which.max(tmp)
          }
        } else if (backfit.method == "dropout") {
          kset <- kset[conv.crit[kset] > backfit.tol]
        }

        if (backfit.method %in% c("extrapolate", "parallel")) {
          if (backfit.method == "extrapolate") {
            proposed.f <- extrapolate.f(flash, old.f, extrapolate.param)
            proposed.f <- update.all.factors(proposed.f)

            old.f <- flash

            if (get.obj(proposed.f) - get.obj(flash) < tol) {
              flash <- update.all.factors(flash)
              extrapolate.param <- decelerate(extrapolate.param)
            } else {
              flash <- proposed.f
              extrapolate.param <- accelerate(extrapolate.param)
            }
          } else if (backfit.method == "parallel") {
            old.f <- flash
            flash <- update.factors.parallel(flash, kset, cl)
          }

          info <- calc.update.info(flash, old.f,
                                   conv.crit.fn, verbose.fns)
          # Since decreases in the objective are possible, use absolute values.
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

        if (is.null(next.tol.target)) {
          next.tol.target <- 10^floor(log10(max(conv.crit)))
        }
        if (max(conv.crit) < next.tol.target) {
          report.backfit.progress(verbose.lvl, next.tol.target)
          next.tol.target <- next.tol.target / 10
        }
      }

      if (backfit.method == "parallel") {
        parallel::stopCluster(cl)
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

      nullchk.tol <- tol * nullchk.reltol
      for (k in nullchk.kset)
        flash <- nullcheck.factor(flash, k, verbose.lvl, nullchk.tol)
      if (nullchk.failed(flash)) {
        if (restart.after.nullchk)
          continue.looping <- TRUE
      } else if (length(nullchk.kset) > 0) {
        report.nullchk.success(verbose.lvl)
      }
    }
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl, is.converged)

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
