# TODO set extrapolate via global options; (also) set verbose in flash obj
flash.add.greedy <- function(flash,
                             Kmax = 1,
                             prior.family = prior.point.normal(),
                             init.fn = init.fn.default,
                             conv.crit.fn = calc.obj.diff,
                             tol = set.default.tol(flash),
                             maxiter = 500,
                             extrapolate = TRUE,
                             verbose.lvl = 1,
                             output.lvl = 3) {
  set.seed(666)

  if (inherits(flash, "flash")) {
    flash <- get.fit(flash)
  }

  priors <- prior.param(prior.family, get.dim(flash))
  ebnm.fns <- rep(priors$ebnm.fn, length.out = Kmax)
  prior.signs <- rep(priors$prior.sign, length.out = Kmax)

  if (missing(tol)) {
    report.tol.setting(verbose.lvl, tol)
  }

  if (uses.R(flash)) {
    if (!missing(extrapolate) && extrapolate) {
      stop("Greedy extrapolation has not yet been implemented for the 'zero' ",
           "variance structure.")
    } else {
      extrapolate <- FALSE
    }
  }

  extrapolate.control <- getOption("extrapolate.control", list())
  extrapolate.param <- set.extrapolate.param(extrapolate.control)

  verbose <- verbose.param(verbose.lvl, get.dim(flash))

  if (verbose.lvl == -1)
    print.tab.delim.table.header(verbose.colnames)

  factors.added <- 0
  greedy.failed <- FALSE
  is.converged <- TRUE

  while (factors.added < Kmax && !greedy.failed) {
    announce.add.factor(verbose.lvl, k = get.next.k(flash))

    factor <- init.factor(flash, init.fn, prior.signs[[factors.added + 1]])
    factor <- set.ebnm.fn(factor, ebnm.fns[[factors.added + 1]])

    announce.factor.opt(verbose.lvl)
    print.table.header(verbose.lvl,
                       verbose$verbose.colnames,
                       verbose$verbose.colwidths,
                       backfit = FALSE)

    iter <- 0
    conv.crit <- Inf
    if (extrapolate) {
      old.f <- factor
      extrapolate.param <- init.beta(extrapolate.param)
    }
    while (conv.crit > tol && iter < maxiter) {
      iter <- iter + 1

      if (extrapolate) {
        proposed.factor <- extrapolate.f(factor, old.f, extrapolate.param)
        proposed.factor <- update.factor(proposed.factor, flash)
      }

      old.f <- factor
      if (!extrapolate) {
        factor <- update.factor(factor, flash)
      } else if (get.obj(proposed.factor) - get.obj(factor) < tol) {
        factor <- update.factor(factor, flash)
        extrapolate.param <- decelerate(extrapolate.param)
      } else {
        factor <- proposed.factor
        extrapolate.param <- accelerate(extrapolate.param)
      }

      obj.diff <- get.obj(factor) - get.obj(old.f)
      if (is.obj.valid(old.f) && obj.diff < 0) {
        report.greedy.obj.decrease(verbose.lvl, obj.diff)
        factor <- old.f
        break
      }

      info <- calc.update.info(factor, old.f, conv.crit.fn, verbose$verbose.fns)
      conv.crit <- get.conv.crit(info)
      print.table.entry(verbose.lvl,
                        verbose$verbose.colwidths,
                        iter,
                        info,
                        get.next.k(flash),
                        backfit = FALSE)
    }

    if (iter == maxiter) {
      is.converged <- FALSE
    }

    if (get.obj(factor) > get.obj(flash) + tol
        || !is.obj.valid(flash, factor)) {
      flash <- add.new.factor.to.flash(factor, flash)
      factors.added <- factors.added + 1
    } else {
      greedy.failed <- TRUE
    }

    report.add.factor.result(verbose.lvl, greedy.failed, get.obj(flash))
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl, is.converged)

  report.completion(verbose.lvl)

  return(flash)
}
