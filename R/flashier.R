flashier <- function(data,
                     S = NULL,
                     flash.init = NULL,
                     var.type = 0,
                     prior.type = "point.normal",
                     greedy.Kmax = 50,
                     backfit = c("none",
                                 "final",
                                 "alternating",
                                 "only"),
                     verbose = 1,
                     ...) {
  data <- set.flash.data(data, S)

  must.be.flash.object(flash.init)
  must.be.valid.integer(var.type, lower = 0, upper = get.dim(data))
  must.be.valid.integer(greedy.Kmax, lower = 0)
  if (!is.character(verbose))
    must.be.valid.integer(verbose, lower = -1, upper = 3)

  workhorse.param <- list()
  ellipsis <- list(...)

  if (is.null(ellipsis$prior.sign) && is.null(ellipsis$ebnm.fn)) {
    workhorse.param <- c(workhorse.param,
                         prior.param(prior.type, get.dim(data), ellipsis))
    ellipsis$ash.param  <- NULL
    ellipsis$ebnm.param <- NULL
  } else if (!missing(prior.type)) {
    stop(paste("If prior.type is specified, then prior.sign and ebnm.fn",
               "cannot be."))
  }

  if (is.null(ellipsis$do.final.backfit)
      && is.null(ellipsis$backfit.after)
      && is.null(ellipsis$backfit.every)) {
    backfit <- match.arg(backfit)
    if (backfit == "only") {
      if (!(missing(greedy.Kmax) || greedy.Kmax == 0))
        stop("Cannot set backfit to only with greedy.Kmax > 0")
      greedy.Kmax <- 0
    }
    workhorse.param <- c(workhorse.param, control.params(backfit))
  } else if (!missing(backfit)) {
    stop(paste("If backfit is specified, then do.final.backfit,",
               "backfit.after, and backfit.every cannot be."))
  }

  if (is.null(ellipsis$verbose.fns)
      && is.null(ellipsis$verbose.colnames)
      && is.null(ellipsis$verbose.colwidths)) {
    workhorse.param <- c(workhorse.param,
                          verbose.params(verbose))
  } else if (is.null(ellipsis$verbose.fns)
             || is.null(ellipsis$verbose.colnames)
             || is.null(ellipsis$verbose.colwidths)) {
    stop(paste("If one of verbose.fns, verbose.colnames, and",
               "verbose.colwidths is specified, then all must be."))
  } else if (!(missing(verbose) || verbose == 3)) {
    stop("Custom verbose output cannot be specified with verbose.lvl < 3")
  } else {
    verbose.lvl <- 3
  }

  return(do.call(flash.workhorse, c(list(Y = get.Y(data),
                                         nonmissing = get.nonmissing(data),
                                         given.tau = get.given.tau(data),
                                         given.tau.dim = get.given.tau.dim(data),
                                         flash.init = flash.init,
                                         est.tau.dim = var.type,
                                         greedy.Kmax = greedy.Kmax),
                                    workhorse.param,
                                    ellipsis)))
}

prior.param <- function(prior.type, data.dim, addl.param) {
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
  prior.sign <- rapply(prior.type, prior.type.to.prior.sign, how = "list")
  ebnm.fn    <- rapply(prior.type, prior.type.to.ebnm.fn, how = "list")
  ebnm.param <- rapply(prior.type, prior.type.to.ebnm.param, how = "list",
                       addl.param = addl.param)

  return(list(prior.sign = prior.sign,
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

prior.type.to.prior.sign <- function(prior.type) {
  return(switch(prior.type,
                nonnegative = 1,
                nonpositive = -1,
                0))
}

prior.type.to.ebnm.fn <- function(prior.type) {
  return(switch(prior.type,
                point.normal =,
                point.laplace =,
                nonzero.mode = ebnm.pn,
                ebnm.ash))
}

prior.type.to.ebnm.param <- function(prior.type, addl.param) {
  ebnm.param <- switch(prior.type,
                       point.normal = list(prior_type = "point_normal"),
                       point.laplace = list(prior_type = "point_laplace"),
                       nonzero.mode = list(prior_type = "point_normal",
                                           fix_mu = FALSE),
                       normal.mixture = list(mixcompdist = "normal"),
                       uniform.mixture = list(mixcompdist = "uniform"),
                       nonnegative = list(mixcompdist = "+uniform"),
                       nonpositive = list(mixcompdist = "-uniform"))

  if (!is.null(ebnm.param[["mixcompdist"]])) {
    # Additional parameters for ashr:
    ebnm.param <- c(ebnm.param, list(method = "shrink", output = "flash_data"))
    ebnm.param <- c(ebnm.param, addl.param$ash.param)
  } else {
    # Additional parameters for ebnm:
    ebnm.param <- c(ebnm.param, addl.param$ebnm.param)
  }

  return(ebnm.param)
}

control.params <- function(backfit) {
  control <- list()
  if (backfit %in% c("final", "only")) {
    control$do.final.backfit <- TRUE
  } else if (backfit == "alternating") {
    control$backfit.after <- 2
    control$backfit.every <- 1
    control$do.final.backfit <- TRUE
  }
  return(control)
}

verbose.params <- function(verbose) {
  verbose.param <- list()
  if (is.character(verbose)) {
    verbose <- unlist(strsplit(verbose, " "))
    verbose.param$verbose.lvl         <- 3
    verbose.param$verbose.fns         <- look.up.verbose.fns(verbose)
    verbose.param$verbose.colnames    <- look.up.verbose.colnames(verbose)
    verbose.param$verbose.colwidths   <- look.up.verbose.colwidths(verbose)
  } else {
    verbose.param$verbose.lvl         <- verbose
    if (verbose > 2) {
      verbose.param$verbose.fns       <- c(calc.obj.diff, calc.max.chg.EF)
      verbose.param$verbose.colnames  <- c("Obj Diff", "Max Chg")
      verbose.param$verbose.colwidths <- c(12, 12)
    }
  }
  return(verbose.param)
}

look.up.verbose.fns <- function(verbose) {
  fns <- lapply(verbose, function(symbol) {
    chars <- unlist(strsplit(symbol, ""))
    if (chars[[1]] == "O")
      return(calc.obj.diff)
    if (chars[[1]] == "L") {
      if (length(chars) > 1) {
        n <- as.integer(chars[[2]])
        return(function(new, old) calc.max.chg.EF(new, old, n))
      } else {
        return(calc.max.chg.EF)
      }
    }
    if (chars[[1]] == "W") {
      if (length(chars) > 1) {
        n <- as.integer(chars[[2]])
        return(function(new, old) which.max.chg.EF(new, old, n))
      } else {
        return(which.max.chg.EF)
      }
    }
    stop("Unrecognized verbose output character.")
  })
  return(fns)
}

look.up.verbose.colnames <- function(verbose) {
  names <- lapply(verbose, function(symbol) {
    chars <- unlist(strsplit(symbol, ""))
    if (chars[[1]] == "O")
      return("Obj Diff")
    if (chars[[1]] == "L") {
      name <- "Max Chg"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    }
    if (chars[[1]] == "W") {
      name <- "Whch"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    }
    stop("Unrecognized verbose output character.")
  })
  return(unlist(names))
}

look.up.verbose.colwidths <- function(verbose) {
  widths <- lapply(verbose, function(symbol) {
    char <- substring(symbol, 1, 1)
    if (char %in% c("O", "L"))
      return(13)
    if (char == "W")
      return(9)
    stop("Unrecognized verbose output character.")
  })
  return(unlist(widths))
}
