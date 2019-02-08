#' Fit an empirical Bayes matrix/tensor factorization
#'
#' TODO: Describe model. Add reference.
#'
#' @param data The observations. Can be a matrix, sparse matrix of class
#'   \code{Matrix}, three-dimensional array, or \code{flash.data} object
#'   obtained from \code{set.flash.data}.
#'
#' @param S The standard errors. Can be a matrix, scalar (if standard errors
#'   are the same for all observations), or vector (if, for example, \code{S}
#'   only varies across columns and is constant within any given row). If
#'   \code{NULL}, all residual variance will be estimated.
#'
#' @param flash.init An initial \code{flash} or \code{flash.fit} object.
#'
#' @param var.type Describes the structure of the estimated residual variance.
#'   Can be \code{NULL}, \code{0}, or a vector. If \code{NULL}, then
#'   \code{S} accounts for all residual variance. Otherwise, a rank-one
#'   variance structure will be estimated (and added to any variance specified
#'   by \code{S}). \code{var.type} then gives the dimensions
#'   along which the residual variance is permitted to vary. For example,
#'   \code{var.type = 1} estimates row-specific residual variances, while
#'   \code{var.type = c(1, 2)} estimates an arbitrary rank-one matrix. If
#'   \code{var.type = 0}, then the residual variance is assumed to be constant
#'   across all observations.
#'
#' @param prior.type Indicates the class of distributions that the priors are
#'   assumed to belong to. Can be a vector of length 1, a vector with one
#'   element for each dimension, or a list of vectors. In the last case,
#'   the first list element gives the class(es) for the first factor, the
#'   second gives the class(es) for the second factor, and so on. The last
#'   list element is then re-used as often as necessary. Options include
#'   \code{"point.normal"}, \code{"point.laplace"}, \code{"normal.mixture"},
#'   \code{"uniform.mixture"}, \code{"nonnegative"}, \code{"nonpositive"}, and
#'   \code{"nonzero.mode"}.
#'
#' @param greedy.Kmax The maximum number of factors to be added. Fixed factors
#'   are not counted.
#'
#' @param backfit Whether and how to backfit. If \code{"final"}, then a single
#'   backfit is performed after as many factors as possible have been added.
#'   If \code{"alternating"}, a backfit will be performed after each
#'   factor is added. Set \code{backfit = "only"} to backfit \code{flash.init}
#'   without adding additional factors.
#'
#' @param ebnm.param Additional parameters to be passed to \code{ebnm::ebnm}.
#'   Used by prior types \code{"point.normal"}, \code{"point.laplace"}, and
#'   \code{"nonzero.mode"}.
#'
#' @param ash.param Additional parameters to be passed to \code{ashr::ash}.
#'   Used by prior types \code{"normal.mixture"}, \code{"uniform.mixture"},
#'   \code{"nonnegative"}, and \code{"nonpositive"}.
#'
#' @param verbose.lvl When and how to display progress updates. Set to
#'   \code{0} for none, \code{1} for updates after a factor is added or a
#'   backfit is completed, \code{2} for additional notifications about the
#'   variational lower bound, and \code{3} for updates after every iteration.
#'   Set to \code{-1} to output a single tab-delimited table of values.
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link{flash.workhorse}}. Not for the faint of heart.
#'
#' @return A \code{flash} object. Contains elements:
#'   \describe{
#'     \item{\code{n.factors}}{The total number of factors in the fitted
#'       model.}
#'     \item{\code{objective}}{The variational lower bound achieved by the
#'       fitted model.}
#'     \item{\code{pve}}{The proportion of variance explained by each factor.}
#'     \item{\code{loadings}}{For each factor, the normalized loadings and a
#'       scaling constant.}
#'     \item{\code{sampler}}{A function that takes a single argument
#'       \code{nsamp} and returns \code{nsamp} samples from the posterior
#'       distribution of the (non-normalized) loadings.}
#'     \item{\code{fit}}{A \code{flash.fit} object.}
#'   }
#'
#' @export
#'
flashier <- function(data,
                     S = NULL,
                     flash.init = NULL,
                     var.type = 0,
                     prior.type = "point.normal",
                     greedy.Kmax = 30,
                     backfit = c("none",
                                 "final",
                                 "alternating",
                                 "only"),
                     ebnm.param = list(),
                     ash.param = list(),
                     verbose.lvl = 1,
                     ...) {
  if (is(data, "flash.data") && !missing(S))
    warning("Data has already been set. Ignoring S.")
  data <- set.flash.data(data, S, var.type = var.type)

  ellipsis <- list(...)

  # When available, use existing flash object settings as defaults.
  if (is(flash.init, "flash"))
    flash.init <- flash.init$fit
  if (!is.null(flash.init) && !is(flash.init, "flash.fit"))
    stop("flash.init must be a flash or flash.fit object.")
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

  # Check arguments.
  must.be.valid.var.type(var.type, get.dim(data))
  must.be.integer(greedy.Kmax, lower = 0)
  must.be.named.list(ebnm.param)
  must.be.named.list(ash.param)
  if (!is.character(verbose.lvl))
    must.be.integer(verbose.lvl, lower = -1, upper = 3)

  workhorse.param <- list()

  # Handle "prior type" parameter.
  if (!is.null(flash.init) && !missing(prior.type) && !is.list(prior.type)) {
    # The last element of ebnm.fn (and ebnm.param) specifies settings for new
    #   factors. If there is an initial flash object, the existing settings
    #   need to be kept, while the last list element is overridden.
    new.prior.param <- prior.param(prior.type,
                                   get.dim(data),
                                   ebnm.param,
                                   ash.param)
    k <- length(flash.init$ebnm.fn)
    flash.init$ebnm.fn[k] <- new.prior.param$ebnm.fn
    flash.init$ebnm.param[k] <- new.prior.param$ebnm.param
    ellipsis$ebnm.fn <- flash.init$ebnm.fn
    ellipsis$ebnm.param <- flash.init$ebnm.param
    ellipsis$prior.sign <- new.prior.param$prior.sign
  } else if (is.null(ellipsis$prior.sign) && is.null(ellipsis$ebnm.fn)) {
    # I can't check that ellipsis$ebnm.param is NULL because I've overloaded
    #   ebnm.param.
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
        stop("Cannot set backfit to \"only\" with greedy.Kmax > 0.")
      greedy.Kmax <- 0
    }
    workhorse.param <- c(workhorse.param, control.param(backfit))
  } else if (!missing(backfit)) {
    stop(paste("If backfit is specified, then final.backfit,",
               "backfit.after, and backfit.every cannot be."))
  }

  # Handle "verbose.lvl" parameter.
  if (is.null(ellipsis$verbose.fns)
      && is.null(ellipsis$verbose.colnames)
      && is.null(ellipsis$verbose.colwidths)) {
    workhorse.param <- c(workhorse.param, verbose.param(verbose.lvl,
                                                        get.dim(data)))
  } else if (is.null(ellipsis$verbose.fns)
             || is.null(ellipsis$verbose.colnames)
             || is.null(ellipsis$verbose.colwidths)) {
    stop(paste("If one of verbose.fns, verbose.colnames, and",
               "verbose.colwidths is specified, then all must be."))
  } else if (!(missing(verbose.lvl) || verbose.lvl %in% c(-1, 3))) {
    stop("Custom verbose output cannot be specified with verbose.lvl < 3")
  } else {
    ellipsis$verbose.lvl <- 3
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

  error.msg <- "Invalid argument to prior.type."

  prior.type <- lapply(prior.type, function(k) {
    if (!is.vector(k))
      stop(error.msg)
    if (length(k) == 1)
      k <- rep(k, data.dim)
    if (length(k) != data.dim)
      stop(error.msg)
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
    # Additional parameters for ashr::ash.
    param <- c(param, list(method = "shrink", output = "flash_data"))
    param <- modifyList(param, ash.param)
  } else {
    # Additional parameters for ebnm::ebnm.
    param <- modifyList(param, ebnm.param)
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
    control$final.backfit <- FALSE
  }
  return(control)
}

verbose.param <- function(verbose, data.dim) {
  param <- list()
  if (is.character(verbose)) {
    verbose <- unlist(strsplit(toupper(verbose), "[ .,/]"))
    verbose <- verbose[verbose != ""]
    param$verbose.lvl         <- 3
    param$verbose.fns         <- look.up.verbose.fns(verbose, data.dim)
    param$verbose.colnames    <- look.up.verbose.colnames(verbose)
    param$verbose.colwidths   <- look.up.verbose.colwidths(verbose)
  } else {
    param$verbose.lvl         <- verbose
    if (verbose > 2) {
      param$verbose.fns       <- c(calc.obj.diff, calc.max.chg.EF)
      param$verbose.colnames  <- c("Obj Diff", "Max Chg")
      param$verbose.colwidths <- c(12, 12)
    } else if (verbose == -1) {
      param$verbose.fns       <- c(get.new.obj, calc.obj.diff, calc.max.chg.EF)
      param$verbose.colnames  <- c("Obj", "Obj.diff", "Max.chg")
      param$verbose.colwidths <- c(14, 12, 12)
    }
  }
  return(param)
}

look.up.verbose.fns <- function(verbose, data.dim) {
  fns <- lapply(verbose, function(symbol) {
    chars <- unlist(strsplit(symbol, ""))
    if (length(chars) > 1) {
      if (length(chars) > 2)
        stop("Unable to parse verbose output string.")
      n <- as.integer(chars[[2]])
      must.be.integer(n, lower = 1, upper = data.dim)
    } else {
      n <- NULL
    }
    if (chars[[1]] == "O") {
      if (!is.null(n))
        warning("Dimension ignored for verbose objective input.")
      return(calc.obj.diff)
    } else if (chars[[1]] == "L") {
      return(function(new, old, k) calc.max.chg.EF(new, old, k, n))
    } else if (chars[[1]] == "W") {
      return(function(new, old, k) which.max.chg.EF(new, old, k, n))
    } else if (chars[[1]] == "S") {
      if (is.null(n))
        stop("Dimension must be specified for verbose sparsity output.")
      return(function(new, old, k) get.sparsity(new, old, k, n))
    } else if (chars[[1]] == "E") {
      return(function(new, old, k) get.exclusion.count(new, old, k, n))
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
    if (chars[[1]] == "S") {
      return(paste("Sparsity", chars[[2]]))
    }
    if (chars[[1]] == "E") {
      return(paste("Excl", chars[[2]]))
    }
    stop("Unrecognized verbose output character.")
  })
  return(unlist(names))
}

look.up.verbose.colwidths <- function(verbose) {
  widths <- lapply(verbose, function(symbol) {
    char <- substring(symbol, 1, 1)
    if (char %in% c("O", "L", "S"))
      return(13)
    if (char %in% c("W", "E"))
      return(9)
    stop("Unrecognized verbose output character.")
  })
  return(unlist(widths))
}
