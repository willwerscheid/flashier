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
                     ebnm.param = NULL,
                     ash.param = NULL,
                     verbose = 1,
                     ...) {
  if (!is(data, "flash.data")) {
    data <- set.flash.data(data, S)
  } else if (!missing(S)) {
    warning("Data has already been set. Ignoring S.")
  }

  ellipsis <- list(...)

  # When available, use existing flash object settings as defaults.
  must.be.flash.object(flash.init)
  if (!is.null(flash.init)) {
    if (missing(var.type))
      var.type <- flash.init$est.tau.dim
    if (missing(prior.type) && is.null(ellipsis$prior.sign))
      ellipsis$prior.sign <- flash.init$dim.signs
    if (missing(prior.type) && is.null(ellipsis$ebnm.fn))
      ellipsis$ebnm.fn <- flash.init$ebnm.fn
    if (missing(prior.type) && is.null(ellipsis$ebnm.param))
      ellipsis$ebnm.param <- flash.init$ebnm.param
  }

  # Check inputs.
  must.be.integer(var.type, lower = 0, upper = get.dim(data))
  must.be.integer(greedy.Kmax, lower = 0)
  must.be.named.list(ebnm.param)
  must.be.named.list(ash.param)
  if (!is.character(verbose))
    must.be.integer(verbose, lower = -1, upper = 3)

  workhorse.param <- list()

  # Handle "prior type" parameter.
  if (is.null(ellipsis$prior.sign) && is.null(ellipsis$ebnm.fn)) {
    workhorse.param <- c(workhorse.param, prior.param(prior.type,
                                                      get.dim(data),
                                                      ebnm.param,
                                                      ash.param))
  } else if (!missing(prior.type)) {
    stop(paste("If prior.type is specified, then prior.sign and ebnm.fn",
               "cannot be."))
  }

  # Handle "backfit" parameter.
  if (is.null(ellipsis$final.backfit)
      && is.null(ellipsis$backfit.after)
      && is.null(ellipsis$backfit.every)) {
    backfit <- match.arg(backfit)
    if (backfit == "only") {
      if (!(missing(greedy.Kmax) || greedy.Kmax == 0))
        stop("Cannot set backfit to only with greedy.Kmax > 0.")
      greedy.Kmax <- 0
    }
    workhorse.param <- c(workhorse.param, control.param(backfit))
  } else if (!missing(backfit)) {
    stop(paste("If backfit is specified, then final.backfit,",
               "backfit.after, and backfit.every cannot be."))
  }

  # Handle "verbose" parameter.
  if (is.null(ellipsis$verbose.fns)
      && is.null(ellipsis$verbose.colnames)
      && is.null(ellipsis$verbose.colwidths)) {
    workhorse.param <- c(workhorse.param, verbose.param(verbose))
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

  return(do.call(flash.workhorse, c(list(data,
                                         flash.init = flash.init,
                                         var.type = var.type,
                                         greedy.Kmax = greedy.Kmax),
                                    workhorse.param,
                                    ellipsis)))
}

prior.param <- function(prior.type, data.dim, ebnm.param, ash.param) {
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
                       ebnm.param = ebnm.param, ash.param = ash.param)

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

prior.type.to.ebnm.param <- function(prior.type, ebnm.param, ash.param) {
  param <- switch(prior.type,
                  point.normal = list(prior_type = "point_normal"),
                  point.laplace = list(prior_type = "point_laplace"),
                  nonzero.mode = list(prior_type = "point_normal",
                                      fix_mu = FALSE),
                  normal.mixture = list(mixcompdist = "normal"),
                  uniform.mixture = list(mixcompdist = "uniform"),
                  nonnegative = list(mixcompdist = "+uniform"),
                  nonpositive = list(mixcompdist = "-uniform"))

  if (!is.null(param[["mixcompdist"]])) {
    # Additional parameters for ashr::ash
    param <- c(param, list(method = "shrink", output = "flash_data"))
    param <- c(param, ash.param)
  } else {
    # Additional parameters for ebnm::ebnm
    param <- c(param, ebnm.param)
  }

  return(param)
}

control.param <- function(backfit) {
  control <- list()
  if (backfit %in% c("final", "only")) {
    control$final.backfit <- TRUE
  } else if (backfit == "alternating") {
    control$backfit.after <- 2
    control$backfit.every <- 1
    control$final.backfit <- TRUE
  }
  return(control)
}

verbose.param <- function(verbose) {
  param <- list()
  if (is.character(verbose)) {
    verbose <- unlist(strsplit(toupper(verbose), "[ .,/]"))
    verbose <- verbose[verbose != ""]
    param$verbose.lvl         <- 3
    param$verbose.fns         <- look.up.verbose.fns(verbose)
    param$verbose.colnames    <- look.up.verbose.colnames(verbose)
    param$verbose.colwidths   <- look.up.verbose.colwidths(verbose)
  } else {
    param$verbose.lvl         <- verbose
    if (verbose > 2) {
      param$verbose.fns       <- c(calc.obj.diff, calc.max.chg.EF)
      param$verbose.colnames  <- c("Obj Diff", "Max Chg")
      param$verbose.colwidths <- c(12, 12)
    }
  }
  return(param)
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
