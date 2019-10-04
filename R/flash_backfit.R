# TODO: need to be able to set prior.family
flash.backfit <- function(flash,
                          kset = NULL,
                          method = c("extrapolate",
                                     "sequential",
                                     "random",
                                     "dropout",
                                     "montaigne",
                                     "parallel"),
                          conv.crit.fn = calc.obj.diff,
                          tol = set.default.tol(flash),
                          maxiter = 500,
                          warmstart = TRUE,
                          verbose.lvl = 1,
                          output.lvl = 3) {
  if (inherits(flash, "flash")) {
    flash <- get.fit(flash)
  }

  if (is.null(kset)) {
    kset <- 1:get.n.factors(flash)
  }

  method <- match.arg(method)

  if (method == "parallel") {
    check.parallel.ok(flash, kset)

    if (missing(conv.crit.fn)) {
      # Since decreases in the ELBO are possible, use absolute values.
      conv.crit.fn <- function(new, old, k) {
        return(abs(calc.obj.diff(new, old, k)))
      }
    }
  }

  if (missing(tol)) {
    report.tol.setting(verbose.lvl, tol)
  }

  # TODO not the smartest way to do this
  flash$warmstart.backfits <- warmstart

  verbose <- verbose.param(verbose.lvl, get.dim(flash))

  conv.crit <- rep(Inf, get.n.factors(flash))
  conv.crit[setdiff(1:get.n.factors(flash), kset)] <- 0

  announce.backfit(verbose.lvl, n.factors = length(kset), tol)
  print.table.header(verbose.lvl,
                     verbose$verbose.colnames,
                     verbose$verbose.colwidths,
                     backfit = TRUE)

  if (method == "parallel") {
    # Remove zero factors and fixed factors.
    kset <- setdiff(kset, which(is.zero(flash)))
    kset <- setdiff(kset, which.k.fixed(flash))

    cl <- parallel::makeCluster(getOption("cl.cores", 2L),
                                type = getOption("cl.type", "PSOCK"),
                                useXDR = FALSE)
  } else if (method == "extrapolate") {
    extrapolate.control <- getOption("extrapolate.control", list())
    extrapolate.param <- set.extrapolate.param(extrapolate.control)
  }

  iter <- 0
  old.obj <- get.obj(flash)
  next.tol.target <- NULL
  if (method == "extrapolate") {
    extrapolate.param <- init.beta(extrapolate.param)
    old.f <- flash
  }
  while (iter < maxiter && max(conv.crit) > tol) {
    is.converged <- TRUE
    iter <- iter + 1

    kset <- get.next.kset(method, kset, conv.crit, tol)

    if (!(method %in% c("parallel", "extrapolate")))  {
      for (k in kset) {
        old.f <- flash
        flash <- update.one.factor(flash, k, iter, verbose.lvl)
        info  <- calc.update.info(flash,
                                  old.f,
                                  conv.crit.fn,
                                  verbose$verbose.fns,
                                  k)
        conv.crit[k] <- get.conv.crit(info)
        print.table.entry(verbose.lvl,
                          verbose$verbose.colwidths,
                          iter,
                          info,
                          k = k,
                          backfit = TRUE)
      }
    } else {
      if (method == "parallel") {
        old.f <- flash
        flash <- update.factors.parallel(flash, kset, cl)
      } else if (method == "extrapolate") {
        proposed.f <- extrapolate.f(flash, old.f, extrapolate.param)
        proposed.f <- update.factors.in.kset(proposed.f, kset)

        old.f <- flash

        if (get.obj(proposed.f) - get.obj(flash) < tol) {
          flash <- update.factors.in.kset(flash, kset)
          extrapolate.param <- decelerate(extrapolate.param)
        } else {
          flash <- proposed.f
          extrapolate.param <- accelerate(extrapolate.param)
        }
      }

      info <- calc.update.info(flash,
                               old.f,
                               conv.crit.fn,
                               verbose$verbose.fns)
      conv.crit <- get.conv.crit(info)
      print.table.entry(verbose.lvl,
                        verbose$verbose.colwidths,
                        iter,
                        info,
                        k = "all",
                        backfit = TRUE)
    }

    if (is.null(next.tol.target) && max(conv.crit) > 0 && max(conv.crit) < Inf) {
      # Set the first target.
      next.tol.target <- 10^floor(log10(max(conv.crit)))
    } else if (max(conv.crit) < next.tol.target) {
      # Report progress and set the next target.
      report.backfit.progress(verbose.lvl, next.tol.target)
      next.tol.target <- next.tol.target / 10
    }
  }

  if (method == "parallel") {
    parallel::stopCluster(cl)
  }

  if (iter == maxiter) {
    is.converged <- FALSE
    report.maxiter.reached(verbose.lvl)
  }

  if (get.obj(flash) > old.obj) {
    report.backfit.complete(verbose.lvl, get.obj(flash))
  }

  # TODO see above
  flash$warmstart.backfits <- NULL

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl, is.converged)

  report.completion(verbose.lvl)
  return(flash)
}

# For parallel backfits to be worthwhile, there must be an efficient way to
#   calculate R2 (that is, more efficient than explicitly forming the n x p
#   matrix of residuals). When the estimated variance is not constant across
#   all entries or if there is any missing data, then (at minimum) a n x k^2
#   or p x k^2 matrix must be created, so it must be true that
#   k^2 << max(n, p) for parallelization to yield any real benefits. But for
#   small k, serial backfits are fast and yield monotonic increases in the
#   ELBO, so parallelization is best avoided.
#
check.parallel.ok <- function(flash, kset) {
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
  if (uses.R(flash)) {
    stop("Parallel backfits cannot be performed with a 'zero' variance ",
         "structure.")
  }
  if (length(which.k.fixed(flash)) > 0) {
    if (length(intersect(which.k.fixed(flash), kset)) > 0) {
      stop("Parallel backfits have not yet been implemented for fixed ",
           "factors.")
    }
    if (is.null(kset)) {
      warning("Parallel backfits have not yet been implemented for fixed ",
              "factors, so they will be removed from kset.")
    }
  }
}

get.next.kset <- function(method, kset, conv.crit, tol) {
  if (method == "random") {
    kset <- sample(kset)
  } else if (method == "montaigne") {
    # Il faut courir au plus pressé.
    last.max <- kset
    kset <- which.max(conv.crit)

    # Don't go after the same factor twice in a row.
    if (identical(last.max, kset)) {
      tmp <- conv.crit
      # Settle for second best.
      tmp[kset] <- 0
      if (max(tmp) > tol) {
        kset <- which.max(tmp)
      }
    }
  } else if (method == "dropout") {
    kset <- kset[conv.crit[kset] > tol]
  }

  return(kset)
}
