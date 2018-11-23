flashier <- function(data,
                     S = NULL,
                     flash.init = NULL,
                     var.type = 0,
                     prior.type = "point.normal",
                     greedy.Kmax = 50,
                     fit.strategy = c("greedy.only",
                                      "backfit.final",
                                      "alternating",
                                      "only.backfit"),
                     verbose.lvl = 1,
                     ...) {
  data <- set.flash.data(data, S)

  must.be.flash.object(flash.init)
  must.be.valid.integer(var.type, lower = 0, upper = get.dim(data))
  must.be.valid.integer(greedy.Kmax, lower = 0)
  must.be.valid.integer(verbose.lvl, lower = -1, upper = 3)

  workhorse.params <- list()

  ellipsis.params <- list(...)
  if (is.null(ellipsis.params$dim.signs)
      && is.null(ellipsis.params$ebnm.fn)
      && is.null(ellipsis.params$ebnm.param)) {
    workhorse.params <- c(workhorse.params,
                          prior.params(prior.type, get.dim(data)))
  } else if (!missing(prior.type)) {
    stop(paste("If prior.type is specified, then dim.signs, ebmn.fn, and",
               "ebnm.param cannot be."))
  }

  if (is.null(ellipsis.params$do.final.backfit)
      && is.null(ellipsis.params$backfit.after)
      && is.null(ellipsis.params$backfit.every)) {
    fit.strategy <- match.arg(fit.strategy)
    if (fit.strategy == "only.backfit") {
      if (!(missing(greedy.Kmax) || greedy.Kmax == 0))
        stop("Cannot set fit.strategy to only.backfit and greedy.Kmax > 0")
      greedy.Kmax <- 0
    }
    workhorse.params <- c(workhorse.params,
                          control.params(fit.strategy))
  } else if (!missing(fit.strategy)) {
    stop(paste("If fit.strategy is specified, then do.final.backfit,",
               "backfit.after, and backfit.every cannot be."))
  }

  if (is.null(ellipsis.params$verbose.fns)
      && is.null(ellipsis.params$verbose.colnames)
      && is.null(ellipsis.params$verbose.colwidths)) {
    workhorse.params <- c(workhorse.params,
                          verbose.params(verbose.lvl))
  } else if (is.null(ellipsis.params$verbose.fns)
             || is.null(ellipsis.params$verbose.colnames)
             || is.null(ellipsis.params$verbose.colwidths)) {
    stop(paste("If one of verbose.fns, verbose.colnames, and",
               "verbose.colwidths is specified, then all must be."))
  }

  return(do.call(flash.workhorse, c(list(Y = get.Y(data),
                                         nonmissing = get.nonmissing(data),
                                         given.tau = get.given.tau(data),
                                         given.tau.dim = get.given.tau.dim(data),
                                         flash.init = flash.init,
                                         est.tau.dim = var.type,
                                         greedy.Kmax = greedy.Kmax,
                                         verbose.lvl = verbose.lvl),
                                    workhorse.params,
                                    ellipsis.params)))
}

prior.params <- function(prior.type, data.dim) {
  if (!is.list(prior.type))
    prior.type <- list(prior.type)

  prior.type <- lapply(prior.type, function(k) {
    if (!is.vector(k))
      stop()
    if (length(k) == 1)
      k <- rep(k, data.dim)
    if (length(k) != data.dim)
      stop()
    return(as.list(k))
  })

  prior.type <- rapply(prior.type, match.prior.type.args, how = "list")
  dim.signs  <- rapply(prior.type, prior.type.to.dim.sign, how = "list")
  ebnm.fn    <- rapply(prior.type, prior.type.to.ebnm.fn, how = "list")
  ebnm.param <- rapply(prior.type, prior.type.to.ebnm.param, how = "list")

  return(list(dim.signs = dim.signs,
              ebnm.fn = ebnm.fn,
              ebnm.param = ebnm.param))
}

match.prior.type.args <- function(prior.type = c("point.normal",
                                                 "point.laplace",
                                                 "normal.mixture",
                                                 "uniform.mixture",
                                                 "nonnegative",
                                                 "nonpositive",
                                                 "nonzero.mode")) {
  return(match.arg(prior.type))
}

prior.type.to.dim.sign <- function(prior.type) {
  return(switch(prior.type,
                nonnegative = 1,
                nonpositive = -1,
                NULL))
}

prior.type.to.ebnm.fn <- function(prior.type) {
  return(switch(prior.type,
                point.normal =,
                point.laplace =,
                nonzero.mode = ebnm.pn,
                ebnm.ash))
}

prior.type.to.ebnm.param <- function(prior.type) {
  ebnm.param <- switch(prior.type,
                       point.normal = list(prior_type = "point_normal"),
                       point.laplace = list(prior_type = "point_laplace"),
                       nonzero.mode = list(prior_type = "point_normal",
                                           fix_mu = FALSE),
                       normal.mixture = list(mixcompdist = "normal"),
                       uniform.mixture = list(mixcompdist = "uniform"),
                       nonnegative = list(mixcompdist = "+uniform"),
                       nonpositive = list(mixcompdist = "-uniform"))

  if (!is.null(ebnm.param[["mixcompdist"]]))
    ebnm.param <- c(ebnm.param, list(method = "shrink", output = "flash_data"))

  return(ebnm.param)
}

control.params <- function(fit.strategy) {
  control <- list()
  if (fit.strategy %in% c("backfit.final", "only.backfit")) {
    control$do.final.backfit <- TRUE
  } else if (fit.strategy == "alternating") {
    control$backfit.after <- 2
    control$backfit.every <- 1
    control$do.final.backfit <- TRUE
  }
  return(control)
}

verbose.params <- function(verbose.lvl) {
  verbose <- list()
  if (verbose.lvl > 2) {
    verbose$verbose.fns       <- c(calc.obj.diff, calc.max.chg.EF)
    verbose$verbose.colnames  <- c("Obj Diff", "Max Chg")
    verbose$verbose.colwidths <- c(12, 12)
  }
  return(verbose)
}
