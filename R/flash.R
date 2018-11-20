# Since updating a factor also updates the matrix of residuals,
#   update.kth.factor needs to be called from within the main loop.
# TODO: use.R can probably be removed; using it almost never seems to improve
#   performance!

flash.workhorse <- function(Y,
                            nonmissing = NULL,
                            F.init = NULL,
                            fix.dim = NULL,
                            fix.idx = NULL,
                            fix.vals = NULL,
                            given.tau = NULL,
                            given.tau.dim = NULL,
                            est.tau.dim = 0,
                            dim.signs = NULL,
                            greedy.ebnm.fn = flashr:::ebnm_pn,
                            greedy.ebnm.param = list(),
                            fix.ebnm.fn = NULL,
                            fix.ebnm.param = NULL,
                            use.fixed.to.est.g = FALSE,
                            use.R = TRUE,
                            greedy.Kmax = 100,
                            backfit.order = c("sequential", "random"),
                            nullchk.fixed = FALSE,
                            backfit.after = NULL,
                            backfit.every = NULL,
                            nullchk.after = backfit.after,
                            nullchk.every = backfit.every,
                            do.final.backfit = FALSE,
                            do.final.nullchk = TRUE,
                            conv.crit.fn = calc.obj.diff,
                            init.maxiter = 100,
                            init.tol = 1e-2,
                            init.verbose = FALSE,
                            greedy.maxiter = 500,
                            greedy.tol = 1e-2,
                            backfit.maxiter = 100,
                            backfit.tol = greedy.tol,
                            final.backfit.maxiter = backfit.maxiter,
                            final.backfit.tol = backfit.tol,
                            verbose.lvl = 1,
                            verbose.fns = NULL,
                            verbose.colnames = NULL,
                            verbose.colwidths = NULL,
                            seed = 666) {
  set.seed(seed)
  backfit.order     <- match.arg(backfit.order)

  announce.flash.init(verbose.lvl)
  flash <- init.flash(Y = Y,
                      nonmissing = nonmissing,
                      F.init = F.init,
                      fix.dim = fix.dim,
                      fix.idx = fix.idx,
                      fix.vals = fix.vals,
                      est.tau.dim = est.tau.dim,
                      given.tau = given.tau,
                      given.tau.dim = given.tau.dim,
                      dim.signs = dim.signs,
                      greedy.ebnm.fn = greedy.ebnm.fn,
                      greedy.ebnm.param = greedy.ebnm.param,
                      fix.ebnm.fn = fix.ebnm.fn,
                      fix.ebnm.param = fix.ebnm.param,
                      use.R = use.R,
                      use.fixed.to.est.g = use.fixed.to.est.g)

  total.factors.added <- 0
  max.factors.to.add  <- greedy.Kmax + get.n.fixed(flash)

  # At least one round of backfitting and nullchecking should be performed
  #   when a non-empty flash object is passed in:
  curr.rnd.factors.added <- get.n.factors(flash)

  # TODO - this kind of parameter handling belongs elsewhere...
  when.to.backfit   <- as.Kset(backfit.after, backfit.every, backfit.maxiter)
  when.to.nullcheck <- as.Kset(nullchk.after, nullchk.every, nullchk.maxiter)

  something.changed <- TRUE

  while (something.changed) {
    something.changed <- FALSE
    flash <- clear.flags(flash)

    greedy.complete <- (total.factors.added >= max.factors.to.add)
    if (!greedy.complete) {
      # TODO: verify works for fixed
      announce.add.factor(verbose.lvl, k = get.next.k(flash))
      announce.factor.init(verbose.lvl)
      factor <- init.factor(flash, init.tol, init.maxiter, init.verbose)

      if (is.fixed(factor)) {
        flash <- add.new.factor.to.flash(factor, flash)
      } else if (greedy.maxiter > 0) {
        announce.factor.opt(verbose.lvl)
        print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths)

        iter <- 0
        conv.crit <- Inf
        while (conv.crit > greedy.tol && iter < greedy.maxiter) {
          iter     <- iter + 1
          old.f    <- factor
          factor   <- update.factor(factor, flash)
          obj.diff <- get.obj(factor) - get.obj(old.f)
          if (is.obj.valid(old.f) && obj.diff < 0) {
            report.greedy.obj.decrease(verbose.lvl, obj.diff)
            factor <- old.f
            break
          }
          info <- calc.update.info(old.f, factor, conv.crit.fn, verbose.fns)
          conv.crit <- get.conv.crit(info)
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info)
        }

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
        total.factors.added    <- total.factors.added + 1
        curr.rnd.factors.added <- curr.rnd.factors.added + 1
      }

      report.add.factor.result(verbose.lvl, failure = greedy.complete)
    }

    if (greedy.complete) {
      do.backfit <- do.final.backfit && (curr.rnd.factors.added > 0)
      do.nullchk <- do.final.nullchk && (curr.rnd.factors.added > 0)
      maxiter    <- final.backfit.maxiter
      tol        <- final.backfit.tol
      curr.rnd.factors.added <- 0
    } else {
      do.backfit <- total.factors.added %in% when.to.backfit
      do.nullchk <- total.factors.added %in% when.to.nullcheck
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
          flash <- update.kth.factor(flash, k, iter, verbose.lvl)
          info  <- calc.update.info(old.f, flash, conv.crit.fn, verbose.fns)
          max.conv.crit <- max(max.conv.crit, get.conv.crit(info))
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info, k)
        }
      }
      # TODO: parallel backfit updates; allow increases by
      #   parallel.monotonicity; might want to do some parallel some seq;
      #   (parallel.iters?) or maybe switch to seq when parallel fails
      #   monotonicity req; in any case flexibility would be good.

      if (get.obj(flash) > old.obj)
        something.changed <- TRUE
    }

    if (do.nullchk) {
      nullchk.kset <- 1:get.n.factors(flash)
      if (!nullchk.fixed)
        nullchk.kset <- setdiff(nullchk.kset, which.k.fixed(flash))

      announce.nullchk(verbose.lvl, n.factors = length(nullchk.kset))

      for (k in nullchk.kset)
        flash <- nullchk.kth.factor(flash, k, verbose.lvl)
      if (nullchk.failed(flash)) {
        something.changed <- TRUE
      } else if (length(nullchk.kset) > 0) {
        report.nullchk.success(verbose.lvl)
      }
    }
  }

  announce.wrapup(verbose.lvl)
  # remove R and then return results

  report.completion(verbose.lvl)

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
