# Since updating a factor also updates the matrix of residuals,
#   update.kth.factor needs to be called from within the main loop.

flash.workhorse <- function(Y,
                            nonmissing = NULL,
                            F.init = NULL,
                            fix.dim = NULL,
                            fix.idx = NULL,
                            fix.vals = NULL,
                            given.tau = NULL,
                            given.tau.dim = NULL,
                            est.tau.dim = 0,
                            greedy.ebnm.fn = flashr:::ebnm_pn,
                            greedy.ebnm.param = list(),
                            fix.ebnm.fn = NULL,
                            fix.ebnm.param = NULL,
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
                            init.maxiter = 100,
                            init.tol = 1e-2,
                            greedy.maxiter = 500,
                            greedy.tol = 1e-2,
                            backfit.maxiter = 100,
                            backfit.tol = greedy.tol,
                            final.backfit.maxiter = backfit.maxiter,
                            final.backfit.tol = backfit.tol,
                            verbose.lvl = 1) {

  announce.flash.init(verbose.lvl)
  flash <- init.flash(Y = Y,
                      nonmissing = nonmissing,
                      F.init = F.init,
                      given.tau = given.tau,
                      given.tau.dim = given.tau.dim,
                      est.tau.dim = est.tau.dim,
                      greedy.ebnm.fn = greedy.ebnm.fn,
                      greedy.ebnm.param = greedy.ebnm.param,
                      fix.ebnm.fn = fix.ebnm.fn,
                      fix.ebnm.param = fix.ebnm.param,
                      use.R = use.R)

  # TODO: set Kmax = greedy.Kmax + factors to add
  total.factors.added <- 0
  max.factors.to.add <- greedy.Kmax
  # Ensure that at least one round of backfitting and nullchecking is
  #   performed when a non-empty flash object is passed in:
  curr.rnd.factors.added <- get.n.factors(flash)

  backfit.order     <- match.arg(backfit.order)
  when.to.backfit   <- as.Kset(backfit.after, backfit.every, backfit.maxiter)
  when.to.nullcheck <- as.Kset(nullchk.after, nullchk.every, nullchk.maxiter)

  something.changed <- TRUE

  while (something.changed) {
    something.changed <- FALSE
    flash <- clear.flags(flash)

    greedy.complete <- (total.factors.added >= max.factors.to.add)
    if (!greedy.complete) {
      # TODO: verify works for fixed
      announce.greedy(verbose.lvl, k = get.n.factors(flash) + 1)
      flash <- add.next.factor(flash,
                               greedy.tol, greedy.maxiter,
                               init.tol, init.maxiter)

      if (greedy.failed(flash)) {
        greedy.complete <- TRUE
      } else {
        something.changed      <- TRUE
        total.factors.added    <- total.factors.added + 1
        curr.rnd.factors.added <- curr.rnd.factors.added + 1
      }

      report.greedy.result(verbose.lvl, failure = greedy.complete)
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

      iter <- 0
      obj.diff <- Inf
      while (iter < maxiter && obj.diff > tol) {
        iter <- iter + 1

        kset <- 1:get.n.factors(flash)
        if (identical(backfit.order, "random"))
          kset <- sample(kset)

        # TODO: other convergence criteria
        old.obj <- get.obj(flash)
        for (k in kset)
          flash <- update.kth.factor(flash)
        obj.diff <- get.obj(flash) - old.obj
      }
      # TODO: parallel backfit updates

      if (iter > 1) {
        something.changed <- TRUE
      }
    }

    if (do.nullchk) {
      nullchk.kset <- 1:get.n.factors(flash)
      if (!nullchk.fixed)
        nullchk.kset <- setdiff(nullchk.kset, which.k.fixed(flash))

      announce.nullchk(verbose.lvl, n.factors = length(nullchk.kset))

      for (k in nullchk.kset)
        flash <- nullchk.kth.factor(flash, k)
      if (nullchk.failed(flash))
        something.changed <- TRUE
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
