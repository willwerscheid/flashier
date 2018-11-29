# TODO: change default backfit.order (and maxiter) after tests have been removed

flash.workhorse <- function(data,
                            flash.init = NULL,
                            var.type = 0,
                            prior.sign = NULL,
                            ebnm.fn = flashr:::ebnm_pn,
                            ebnm.param = list(),
                            greedy.Kmax = 100,
                            backfit.order = c("sequential",
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
                            conv.crit.fn = calc.obj.diff,
                            verbose.lvl = 1,
                            verbose.fns = NULL,
                            verbose.colnames = NULL,
                            verbose.colwidths = NULL,
                            output.lvl = 2,
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
                            backfit.maxiter = 100,
                            backfit.tol = greedy.tol,
                            inner.backfit.maxiter = backfit.maxiter,
                            inner.backfit.tol = backfit.tol,
                            seed = 666,
                            use.R = FALSE) {
  set.seed(seed)
  backfit.order <- match.arg(backfit.order)
  if (force.use.R(data, var.type)) {
    if (!missing(use.R) && !use.R)
      stop("R must be used with the requested variance structure.")
    use.R <- TRUE
  }

  announce.flash.init(verbose.lvl)
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
                      use.R = use.R)

  if (is.null(greedy.tol)) {
    greedy.tol <- sqrt(.Machine$double.eps) * prod(get.dims(flash))
    report.tol.setting(verbose.lvl, greedy.tol)
  }

  total.factors.added <- get.n.factors(flash)
  max.factors.to.add  <- greedy.Kmax + get.n.fixed(flash)
  when.to.backfit <- as.Kset(backfit.after, backfit.every, max.factors.to.add)
  when.to.nullchk <- as.Kset(nullchk.after, nullchk.every, max.factors.to.add)

  # At least one round of backfitting and nullchecking should be performed
  #   when a non-empty flash object is passed in.
  curr.rnd.factors.added <- get.n.factors(flash)

  if (verbose.lvl == -1)
    print.tab.delim.table.header(verbose.colnames)

  something.changed <- TRUE
  is.converged      <- TRUE

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
        print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                           backfit = FALSE)

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
          print.table.entry(verbose.lvl, verbose.colwidths, iter, info,
                            get.next.k(flash), backfit = FALSE)
        }

        if (iter == greedy.maxiter)
          is.converged <- FALSE

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

      report.add.factor.result(verbose.lvl, greedy.complete, get.obj(flash))
    }

    if (greedy.complete) {
      do.backfit  <- final.backfit && (curr.rnd.factors.added > 0)
      do.nullchk  <- final.nullchk && (curr.rnd.factors.added > 0)
      maxiter     <- backfit.maxiter
      tol         <- backfit.tol
      curr.rnd.factors.added <- 0
    } else {
      do.backfit  <- total.factors.added %in% when.to.backfit
      do.nullchk  <- total.factors.added %in% when.to.nullchk
      maxiter     <- inner.backfit.maxiter
      tol         <- inner.backfit.tol
    }

    if (do.backfit) {
      announce.backfit(verbose.lvl, n.factors = get.n.factors(flash))
      print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                         backfit = TRUE)

      kset <- 1:get.n.factors(flash)

      iter <- 0
      old.obj <- get.obj(flash)
      conv.crit <- rep(Inf, get.n.factors(flash))
      while (iter < maxiter && max(conv.crit) > tol) {
        is.converged <- TRUE
        iter <- iter + 1

        if (backfit.order == "random") {
          kset <- sample(kset)
        } else if (backfit.order == "montaigne") {
          # Il faut courir au plus pressé.
          kset <- which.max(conv.crit)
        }

        if (backfit.order == "parallel") {
          old.f <- flash
          flash <- update.factors.parallel(flash, kset)
          info  <- calc.update.info(flash, old.f,
                                    conv.crit.fn, verbose.fns)
          max.conv.crit <- abs(get.conv.crit(info))
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
                              k, backfit = TRUE)
          }
        }
      }

      if (iter == maxiter)
        is.converged <- FALSE

      if (get.obj(flash) > old.obj) {
        report.backfit.complete(verbose.lvl, get.obj(flash))
        something.changed <- TRUE
      }
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
  return(!(is.tau.simple(tmp) || is.var.type.kronecker(tmp)))
}

as.Kset <- function(after, every, maxiter) {
  if (is.null(every))
    return(after)
  max.after <- max(c(0, after))
  n.every <- floor((maxiter - max.after) / every)
  every.iters <- max.after + every * (1:n.every)
  return(c(after, every.iters))
}
