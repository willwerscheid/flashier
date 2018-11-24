# Since updating a factor also updates the matrix of residuals,
#   update.kth.factor needs to be called from within the main loop.
# TODO: think more about how to set default tols.
#   Maybe sqrt(machine.eps) * n * p??

flash.workhorse <- function(Y,
                            nonmissing = NULL,
                            given.tau = NULL,
                            given.tau.dim = NULL,
                            flash.init = NULL,
                            est.tau.dim = 0,
                            prior.sign = NULL,
                            candidate.factors = NULL,
                            ebnm.fn = flashr:::ebnm_pn,
                            ebnm.param = list(),
                            greedy.Kmax = 100,
                            backfit.after = NULL,
                            backfit.every = NULL,
                            do.final.backfit = FALSE,
                            backfit.order = c("sequential", "random"),
                            warmstart.backfits = TRUE,
                            nullchk.after = NULL,
                            nullchk.every = NULL,
                            do.final.nullchk = TRUE,
                            conv.crit.fn = calc.obj.diff,
                            verbose.lvl = 1,
                            verbose.fns = NULL,
                            verbose.colnames = NULL,
                            verbose.colwidths = NULL,
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
                            greedy.tol = 1e-2,
                            backfit.maxiter = 100,
                            backfit.tol = greedy.tol,
                            final.backfit.maxiter = backfit.maxiter,
                            final.backfit.tol = backfit.tol,
                            seed = 666,
                            use.R = FALSE) {
  set.seed(seed)
  backfit.order <- match.arg(backfit.order)
  when.to.backfit <- as.Kset(backfit.after, backfit.every, backfit.maxiter)
  when.to.nullchk <- as.Kset(nullchk.after, nullchk.every, nullchk.maxiter)
  if (force.use.R(given.tau)) {
    if (!missing(use.R))
      stop("R must be used if S is a matrix")
    use.R <- TRUE
  }

  announce.flash.init(verbose.lvl)
  if (is.null(flash.init)) {
    flash <- init.new.flash(Y = Y,
                            nonmissing = nonmissing,
                            given.tau = given.tau,
                            given.tau.dim = given.tau.dim,
                            EF.init = EF.init,
                            est.tau.dim = est.tau.dim,
                            dim.signs = prior.sign,
                            ebnm.fn = ebnm.fn,
                            ebnm.param = ebnm.param,
                            warmstart.backfits = warmstart.backfits,
                            fix.dim = fix.dim,
                            fix.idx = fix.idx,
                            fix.vals = fix.vals,
                            use.fixed.to.est.g = use.fixed.to.est.g,
                            use.R = use.R)
  } else {
    # TODO change non-missing arguments in flash.init while arg checking
    flash <- flash.init
  }

  total.factors.added <- get.n.factors(flash)
  max.factors.to.add  <- greedy.Kmax + get.n.fixed(flash)

  # At least one round of backfitting and nullchecking should be performed
  #   when a non-empty flash object is passed in:
  curr.rnd.factors.added <- get.n.factors(flash)

  something.changed <- TRUE

  while (something.changed) {
    something.changed <- FALSE
    flash <- clear.flags(flash)

    greedy.complete <- (total.factors.added >= max.factors.to.add)
    if (!greedy.complete) {
      announce.add.factor(verbose.lvl, k = get.next.k(flash))

      announce.factor.init(verbose.lvl)
      factor <- init.factor(flash, init.fn, init.tol, init.maxiter)

      if (is.fixed(factor)) {
        flash <- add.new.factor.to.flash(factor, flash)
      } else if (greedy.maxiter > 0) {
        announce.factor.opt(verbose.lvl)
        print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths)

        iter <- 0
        conv.crit <- Inf
        while (conv.crit > greedy.tol && iter < greedy.maxiter) {
          iter <- iter + 1

          old.f    <- factor
          factor   <- update.factor(factor, flash)
          obj.diff <- get.obj(factor) - get.obj(old.f)
          if (is.obj.valid(old.f) && obj.diff < 0) {
            report.greedy.obj.decrease(verbose.lvl, obj.diff)
            factor <- old.f
            break
          }
          info <- calc.update.info(factor, old.f, conv.crit.fn, verbose.fns)
          conv.crit <- get.conv.crit(info)
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info)
        }

        # TODO if "same signs" factor, do a "delayed add"; that is, save the
        #   factor in some candidate.factor variable (which can actually
        #   be a list of multiple factors), set something.changed
        #   to TRUE, and compare all factors when we get to the last one

        # if (is.candidate(factor))
        #   candidate.factors <- c(candidate.factors, list(factor))
        # else
        #   if (!is.null(candidate.factors))
        #     factor <- select.best.factor(factor, candidate.factors)
        #     candidate.factors <- NULL

        # this whole thing is part of the else:
        if (get.obj(factor) > get.obj(flash) || !is.obj.valid(flash, factor)) {
          flash <- add.new.factor.to.flash(factor, flash)
        } else {
          flash <- set.greedy.fail.flag(flash)
        }
      }

      if (greedy.failed(flash)) {
        greedy.complete <- TRUE
      } else {
        something.changed      <- TRUE
        # if (!is.candidate):
        total.factors.added    <- total.factors.added + 1
        curr.rnd.factors.added <- curr.rnd.factors.added + 1
      }

      report.add.factor.result(verbose.lvl, failure = greedy.complete)
    }

    # if !is.null(candidate.factors) FALSE, FALSE
    if (greedy.complete) {
      do.backfit <- do.final.backfit && (curr.rnd.factors.added > 0)
      do.nullchk <- do.final.nullchk && (curr.rnd.factors.added > 0)
      maxiter    <- final.backfit.maxiter
      tol        <- final.backfit.tol
      curr.rnd.factors.added <- 0
    } else {
      do.backfit <- total.factors.added %in% when.to.backfit
      do.nullchk <- total.factors.added %in% when.to.nullchk
      maxiter    <- backfit.maxiter
      tol        <- backfit.tol
    }

    if (do.backfit) {
      announce.backfit(verbose.lvl, n.factors = get.n.factors(flash))
      print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                         backfit = TRUE)

      iter <- 0
      max.conv.crit <- Inf
      old.obj <- get.obj(flash)
      while (iter < maxiter && max.conv.crit > tol) {
        iter <- iter + 1
        max.conv.crit <- 0

        kset <- 1:get.n.factors(flash)
        if (identical(backfit.order, "random"))
          kset <- sample(kset)
        for (k in kset) {
          old.f <- flash
          flash <- update.existing.factor(flash, k, iter, verbose.lvl)
          info  <- calc.update.info(flash, old.f, conv.crit.fn, verbose.fns)
          max.conv.crit <- max(max.conv.crit, get.conv.crit(info))
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info, k)
        }
      }
      # TODO: parallel backfit updates; allow increases by some
      #   parallel.monotonicity parameter; might want to do some parallel,
      #   some seq, or maybe switch to seq when parallel fails
      #   monotonicity tol; in any case flexibility would be good.

      if (get.obj(flash) > old.obj)
        something.changed <- TRUE
    }

    if (do.nullchk) {
      nullchk.kset <- 1:get.n.factors(flash)
      if (!nullchk.fixed.factors)
        nullchk.kset <- setdiff(nullchk.kset, which.k.fixed(flash))

      announce.nullchk(verbose.lvl, n.factors = length(nullchk.kset))

      for (k in nullchk.kset)
        flash <- nullcheck.factor(flash, k, verbose.lvl)
      if (nullchk.failed(flash)) {
        something.changed <- TRUE
      } else if (length(nullchk.kset) > 0) {
        report.nullchk.success(verbose.lvl)
      }
    }
  }

  announce.wrapup(verbose.lvl)
  # TODO: remove R and then return results; propagate names here;
  #   include sampler

  report.completion(verbose.lvl)

  class(flash) = "flash"

  return(flash)
}

as.Kset <- function(after, every, maxiter) {
  if (is.null(every))
    return(after)
  max.after <- max(c(0, after))
  n.every <- floor((maxiter - max.after) / every)
  every.iters <- max.after + every * (1:n.every)
  return(c(after, every.iters))
}
