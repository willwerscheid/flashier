#' Backfit a flash object
#'
#' Backfits existing flash factor/loadings pairs. Whereas a "greedy" fit optimizes
#'   each newly added factor/loadings pair in one go without returning to optimize
#'   previously added pairs, a "backfit" updates all existing pairs in a cyclical
#'   fashion. See \code{\link{flash}} for examples of usage.
#'
#' @inheritParams flash
#'
#' @inheritParams flash_greedy
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers specifying which factors to backfit.
#'   If \code{kset = NULL}, then all existing factors will be backfitted.
#'
#' @param maxiter The maximum number of backfitting iterations. An "iteration"
#'   is defined such that all factors in \code{kset} get updated at each
#'   iteration.
#'
#' @param tol The convergence tolerance parameter. After each update, the fit
#'   is compared to the fit from before the update using a convergence
#'   criterion function (by default, the difference in ELBO, but the criterion
#'   can be changed via \code{\link{flash_set_conv_crit}}).
#'   The backfit is considered to have "converged" when the value of the
#'   convergence criterion function over successive updates to
#'   \emph{all} factor/loadings pairs is less than or equal to \code{tol}. If,
#'   for example, factor/loadings pairs \eqn{1, \ldots, K} are being
#'   sequentially backfitted, then fits are compared before and
#'   after the update to factor/loadings 1, before and after the update to
#'   factor/loadings 2, and so on through factor/loadings \eqn{K},
#'   and backfitting only terminates when the convergence criterion function
#'   returns a value less
#'   than or equal to \code{tol} for all \eqn{K} updates. Note that
#'   specifying \code{tol} here will override any value set by
#'   \code{flash_set_conv_crit}; to use the "global" tolerance parameter,
#'   \code{tol} must be left unspecified (\code{NULL}).
#'   If \code{tol = NULL} and a global tolerance
#'   parameter has not been set, then the default
#'   tolerance used is \eqn{np\sqrt{\epsilon}}, where \eqn{n} is the
#'   number of rows in the dataset, \eqn{p} is the number of columns, and
#'   \eqn{\epsilon} is equal to \code{\link{.Machine}$double.eps}.
#'
#' @importFrom parallel makeCluster stopCluster
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, backfitted
#'   as specified.
#'
#' @export
#'
flash_backfit <- function(flash,
                          kset = NULL,
                          extrapolate = TRUE,
                          warmstart = TRUE,
                          maxiter = 500,
                          tol = NULL,
                          verbose = NULL) {
  flash <- get.fit(flash)

  tol <- handle.tol.param(tol, flash)
  verbose.lvl <- handle.verbose.param(verbose, flash)

  if (is.timed.out(flash)) {
    report.timeout.no.backfit(verbose.lvl)
    verbose.lvl <- 0
  }

  if (is.null(kset)) {
    if (get.n.factors(flash) > 0) {
      kset <- 1:get.n.factors(flash)
    } else {
      announce.no.backfit(verbose.lvl)
      verbose.lvl <- 0
    }
  } else {
    must.be.valid.kset(flash, kset)
  }

  must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)

  # Removed methods "dropout", "random", and "parallel".
  if (extrapolate) {
    method <- "extrapolate"
  } else {
    method <- "sequential"
  }

  if (method == "parallel") {
    check.parallel.ok(flash, kset)

    # TODO: address this
    if (missing(conv.crit.fn)) {
      # Since decreases in the ELBO are possible, use absolute values.
      conv.crit.fn <- function(new, old, k) {
        return(abs(calc.obj.diff(new, old, k)))
      }
    }
  }

  flash <- set.warmstart(flash, warmstart)

  verbose.fns <- get.verbose.fns(flash)
  verbose.colnames <- get.verbose.colnames(flash)
  verbose.colwidths <- get.verbose.colwidths(flash)

  conv.crit.fn <- get.conv.crit.fn(flash)

  conv.crit <- rep(Inf, get.n.factors(flash))
  conv.crit[setdiff(1:get.n.factors(flash), kset)] <- 0

  announce.backfit(verbose.lvl, n.factors = length(kset), tol)
  print_table.header(verbose.lvl,
                     verbose.colnames,
                     verbose.colwidths,
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
  while (iter < maxiter && max(conv.crit) > tol && !is.timed.out(flash)) {
    iter <- iter + 1

    kset <- get.next.kset(method, kset, conv.crit, tol)

    if (!(method %in% c("parallel", "extrapolate")))  {
      for (k in kset) {
        old.f <- flash
        flash <- update.one.factor(flash, k, iter, verbose.lvl)
        info  <- calc.update.info(flash,
                                  old.f,
                                  conv.crit.fn,
                                  verbose.fns,
                                  k)
        conv.crit[k] <- get.conv.crit(info)
        print_table.entry(verbose.lvl,
                          verbose.colwidths,
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
                               verbose.fns)
      conv.crit <- get.conv.crit(info)
      print_table.entry(verbose.lvl,
                        verbose.colwidths,
                        iter,
                        info,
                        k = "all",
                        backfit = TRUE)
    }

    if (is.null(next.tol.target) && max(conv.crit) > 0 && max(conv.crit) < Inf) {
      # Set the first target.
      next.tol.target <- 10^floor(log10(max(conv.crit)))
    } else if (!is.null(next.tol.target) && max(conv.crit) < next.tol.target) {
      # Report progress and set the next target.
      report.backfit.progress(verbose.lvl, next.tol.target)
      next.tol.target <- next.tol.target / 10
    }
  }

  if (iter > 0 && is.timed.out(flash)) {
    t.diff <- Sys.time() - get.timeout.set.time(flash)
    report.timeout.reached(verbose.lvl, t.diff)
    flash <- set.timeout.reached.flag(flash)
  }

  if (iter == maxiter) {
    report.maxiter.reached(verbose.lvl)
    flash <- set.max.backfit.iter.reached.flag(flash)
  } else {
    flash <- clear.max.backfit.iter.reached.flag(flash)
  }

  if (method == "parallel") {
    parallel::stopCluster(cl)
  }

  if (get.obj(flash) > old.obj) {
    report.backfit.complete(verbose.lvl, get.obj(flash))
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl = 3L)

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
  if (any_missing(flash)) {
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
