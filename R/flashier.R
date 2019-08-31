#' Empirical Bayes matrix/tensor factorization
#'
#' TODO: Describe model. Add reference.
#'
#' @param data The observations. Can be a matrix, sparse matrix of class
#'   \code{Matrix}, three-dimensional array, or \code{flash.data} object
#'   obtained from \code{\link{set.flash.data}}. Can also be a low-rank matrix
#'   representation as returned by, for example, \code{\link{svd}},
#'   \code{\link[irlba]{irlba}}, \code{\link[rsvd]{rsvd}}, or
#'   \code{\link[softImpute]{softImpute}}. Can be \code{NULL} if
#'   \code{init} is used.
#'
#' @param S The standard errors. Can be a matrix, scalar (if standard errors
#'   are the same for all observations), or vector (if, for example, \code{S}
#'   only varies across columns and is constant within any given row). If
#'   \code{NULL}, all residual variance will be estimated.
#'
#' @param prior.family Indicates the family of distributions that the priors on
#'   the loadings are
#'   assumed to belong to. Can be a list of length 1 or length \eqn{N}, where
#'   \eqn{N} is the number of modes (\eqn{N = 2} for matrices; \eqn{N = 3} for
#'   tensors). Each list element must be a prior family defined by one of the
#'   convenience functions \code{\link{prior.normal}},
#'   \code{\link{prior.point.normal}},
#'   \code{\link{prior.point.laplace}}, \code{\link{prior.nonzero.mode}},
#'   \code{\link{prior.scale.normal.mix}}, \code{\link{prior.unimodal}},
#'   \code{\link{prior.symmetric.unimodal}},
#'   \code{\link{prior.nonnegative}}, or \code{\link{prior.nonpositive}},
#'   or a custom prior type of a similar form (see \code{\link{prior.normal}}
#'   for details).
#'   For example, the default \code{prior.family = prior.point.normal()} fits a
#'   (different) point-normal prior for each factor and each mode, while
#'   \code{prior.family = c(prior.nonnegative(), prior.scale.normal.mix())}
#'   fits a unimodal distribution with mode zero and nonnegative support to
#'   each set of row loadings and a scale mixture of normals with mean zero to
#'   each set of column loadings.
#'
#'   \code{prior.family} can also be a list of lists, in which case the first
#'   list specifies the family or families for the first factor, the second
#'   specifies the family or families for the second factor, and so on. The
#'   last list element is then re-used as often as necessary.
#'   For example, \code{prior.family = list(prior.nonzero.mode(),
#'   prior.scale.normal.mix())} will fit point-normal priors with nonzero
#'   means for the first factor and scale mixtures of normals for every
#'   subsequent factor.
#'
#' @param var.type Describes the structure of the estimated residual variance.
#'   Can be \code{NULL}, \code{0}, or a vector. If \code{NULL}, then
#'   \code{S} accounts for all residual variance. Otherwise, a rank-one
#'   variance structure will be estimated (and added to any variance specified
#'   by \code{S}). \code{var.type} then gives the modes along which the
#'   residual variance is permitted to vary. For example, \code{var.type = 1}
#'   estimates row-specific residual variances, while \code{var.type = c(1, 2)}
#'   optimizes over all rank-one matrices. If \code{var.type = 0}, then the
#'   residual variance is assumed to be constant across all observations.
#'
#' @param fit Fitting method. When \code{fit = "greedy"}, \code{flashier}
#'   adds as many as \code{greedy.Kmax} factors, optimizing each newly added
#'   factor in one go without returning to optimize previously added factors.
#'   When \code{fit = "full"}, \code{flashier} will perform a final "backfit"
#'   where all factors are cyclically updated until convergence. Set
#'   \code{fit = "backfit.only"} to backfit \code{init} without greedily
#'   adding factors.
#'
#' @param init An initial \code{flash} or \code{flash.fit} object.
#'
#' @param fixed.factors Adds factors with fixed loadings. Options
#'   include mean factors (where all row or column loadings are fixed at 1),
#'   factors with known sparsity patterns, and factors with arbitrarily fixed
#'   elements. See \code{\link{fixed.ones}}, \code{\link{fixed.sparse}},
#'   and \code{\link{fixed.factors}} for usage. Multiple types of fixed factors
#'   can be added by concatenating via \code{c()}. For example,
#'   \code{fixed.factors = c(fixed.ones(n = 1), fixed.sparse(n = 1,
#'   nz.idx = 1:10)} will add one mean factor and one sparse factor.
#'
#' @param greedy.Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by \code{flashier}, since
#'   factors are only added as long as they increase the variational lower
#'   bound on the log likelihood for the model. Fixed factors are not counted
#'   towards this limit.
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
#'     \item{\code{pve}}{The proportion of variance explained by each factor.}
#'     \item{\code{loadings.scale, loadings.pm}}{Posterior means for loadings.
#'       Since the model is not identifiable,
#'       each column of loadings is \eqn{L2}-normalized. The
#'       normalization constant is given by \code{loadings.scale}. Thus, for
#'       matrices, fitted values can be calculated as \code{f$loadings.pm[[1]]
#'       \%*\% diag(f$loadings.scale) \%*\% t(f$loadings.pm[[2]])} (or, more
#'       simply, as \code{fitted(f)}).}
#'     \item{\code{loadings.psd}}{Posterior standard deviations for loadings.}
#'     \item{\code{loadings.lfsr}}{Local false sign rates for loadings.}
#'     \item{\code{residuals.sd}}{Estimated residual standard deviations (these
#'       include any variance component given as an argument to \code{S}).}
#'     \item{\code{fitted.g}}{The fitted priors for each mode and factor.}
#'     \item{\code{elbo}}{The variational lower bound achieved by the
#'       fitted model.}
#'     \item{\code{convergence.status}}{A character string indicating whether
#'       the fitting algorithm has converged.}
#'     \item{\code{sampler}}{A function that takes a single argument
#'       \code{nsamp} and returns \code{nsamp} samples from the posterior
#'       distribution of the (non-normalized) loadings.}
#'     \item{\code{flash.fit}}{A \code{flash.fit} object. Used by
#'       \code{flashier} when fitting is not performed all at once, but
#'       incrementally via repeated calls to \code{flashier} (with the
#'       intermediate \code{flash} or \code{flash.fit} objects given as
#'       arguments to \code{init}).}
#'   }
#'
#' @examples
#' # TODO
#'
#' @export
#'
flashier <- function(data = NULL,
                     S = NULL,
                     prior.family = prior.point.normal(),
                     var.type = 0L,
                     fit = c("greedy",
                             "full",
                             "backfit.only"),
                     init = NULL,
                     fixed.factors = NULL,
                     greedy.Kmax = 50L,
                     verbose.lvl = 1L,
                     ...) {
  if (inherits(data, "flash.data") && !missing(S))
    warning("Data has already been set. Ignoring S.")

  ellipsis <- list(...)

  if (!missing(prior.family)
      && (!is.null(ellipsis$prior.sign)
          || !is.null(ellipsis$ebnm.fn)))
    stop(paste("If prior.family is specified, then prior.sign and ebnm.fn ",
               "cannot be."))

  # When available, use existing flash object settings as defaults.
  if (inherits(init, "flash"))
    init <- get.fit(init)
  if (!is.null(init)) {
    if (!inherits(init, "flash.fit"))
      stop("init must be a flash or flash.fit object.")
    if (missing(var.type)) {
      var.type <- get.est.tau.dim(init)
    } else if (!identical(var.type, get.est.tau.dim(init))) {
      init <- clear.bypass.init.flag(init)
      if (is.null(data))
        stop("When changing var.type, data cannot be NULL.")
    }
    if (missing(prior.family)) {
      if (is.null(ellipsis$prior.sign))
        ellipsis$prior.sign <- init$dim.signs
      if (is.null(ellipsis$ebnm.fn))
        ellipsis$ebnm.fn <- init$ebnm.fn
    }
  }

  # Bypass set.flash.data if init has the needed fields.
  if (bypass.init(init)) {
    if (!is.null(data) || !is.null(S))
      warning("Flash object does not need to be re-initialized. Ignoring",
              " data (and/or S).")
    data <- NULL
    dims <- get.dims(init)
  } else {
    data <- set.flash.data(data, S, var.type = var.type)
    dims <- get.dims(data)
  }
  data.dim <- length(dims)

  # Check arguments.
  must.be.valid.var.type(var.type, data.dim)
  must.be.integer(greedy.Kmax, lower = 0)
  if (!is.null(fixed.factors))
    must.be.list.of.named.lists(fixed.factors, c("dim", "idx", "vals"))
  if (!is.character(verbose.lvl))
    must.be.integer(verbose.lvl, lower = -1, upper = 3)

  workhorse.param <- list()

  # Handle "prior family" parameter.
  if (is.null(init)) {
    workhorse.param <- c(workhorse.param, prior.param(prior.family, data.dim))
  } else if (!missing(prior.family)) {
    # The last element of ebnm.fn specifies settings for new factors. If there
    #   is an initial flash object, the existing settings need to be kept,
    #   while the last list element is overridden.
    k <- length(init$ebnm.fn)
    init$dim.signs <- init$dim.signs[-k]
    init$ebnm.fn <- init$ebnm.fn[-k]
    new.prior.param <- prior.param(prior.family, data.dim)
    ellipsis$prior.sign <- c(init$dim.signs, new.prior.param$prior.sign)
    ellipsis$ebnm.fn <- c(init$ebnm.fn, new.prior.param$ebnm.fn)
  }

  # Handle "fixed factors" parameter.
  if (length(fixed.factors) > 0) {
    fixed.factors <- lapply(fixed.factors, function(f) {
      # Handle fixed factor modes.
      dim <- f$dim
      must.be.integer(dim, lower = 1, upper = data.dim, allow.null = FALSE)
      # Handle fixed factor indices (the default is to fix the entire factor).
      idx <- 1:dims[dim]
      if (!is.null(f$idx) && length(f$idx) > 0) {
        if (!all(f$idx %in% c(-idx, idx)))
          stop("Invalid fixed factor indices.")
        idx <- idx[f$idx]
      }
      # Handle fixed factor values.
      vals <- f$vals
      if (length(vals) == 1) {
        vals <- rep(vals, length(idx))
      } else if (length(vals) != length(idx)) {
        stop("The lengths of the fixed factor indices and values do not",
             " match.")
      }
      return(list(dim = dim, idx = idx, vals = vals))
    })

    # If init is used, previously fixed factors need to be included.
    #   These fields aren't always populated, so an offset needs to be
    #   used when setting the new fixed factors.
    ellipsis$fix.dim <- list()
    ellipsis$fix.idx <- list()
    ellipsis$fix.vals <- list()
    if (!is.null(get.fix.dim(init))) {
      ellipsis$fix.dim  <- get.fix.dim(init)
      ellipsis$fix.idx  <- get.fix.idx(init)
      ellipsis$fix.vals <- get.fix.vals(init)
    }
    kset <- get.n.factors(init) + 1:length(fixed.factors)
    ellipsis$fix.dim[kset]  <- lapply(fixed.factors, `[[`, "dim")
    ellipsis$fix.idx[kset]  <- lapply(fixed.factors, `[[`, "idx")
    ellipsis$fix.vals[kset] <- lapply(fixed.factors, `[[`, "vals")
  }

  # Handle "fit" parameter.
  if (is.null(ellipsis$final.backfit)) {
    fit <- match.arg(fit)
    if (fit == "greedy") {
      ellipsis$final.backfit <- FALSE
    } else if (fit == "backfit.only") {
      if (is.null(init) && is.null(fixed.factors) && is.null(ellipsis$fix.dim)
          && is.null(ellipsis$EF.init))
        stop("There's nothing to backfit. Did you mean to set fit = \"full\"?")
      if (!(missing(greedy.Kmax) || greedy.Kmax == 0))
        stop("Cannot set fit to \"backfit.only\" with greedy.Kmax > 0.")
      greedy.Kmax <- 0
      ellipsis$final.backfit <- TRUE
    } else if (fit == "full") {
      ellipsis$final.backfit <- TRUE
    }
  } else if (!missing(fit)) {
    stop("If fit is specified, then final.backfit cannot be.")
  }

  # Handle "verbose.lvl" parameter.
  if (is.null(ellipsis$verbose.fns)
      && is.null(ellipsis$verbose.colnames)
      && is.null(ellipsis$verbose.colwidths)) {
    workhorse.param <- c(workhorse.param, verbose.param(verbose.lvl, data.dim))
  } else if (is.null(ellipsis$verbose.fns)
             || is.null(ellipsis$verbose.colnames)
             || is.null(ellipsis$verbose.colwidths)) {
    stop(paste("If one of verbose.fns, verbose.colnames, and",
               "verbose.colwidths is specified, then all must be."))
  } else if (missing(verbose.lvl)) {
    ellipsis$verbose.lvl <- 3
  } else if (!(verbose.lvl %in% c(-1, 3))) {
    stop("If custom verbose output is specified, then verbose.lvl must be -1",
         " or 3.")
  } else {
    ellipsis$verbose.lvl <- verbose.lvl
  }

  return(do.call(flash.workhorse, c(list(data,
                                         init = init,
                                         var.type = var.type,
                                         greedy.Kmax = greedy.Kmax),
                                    workhorse.param,
                                    ellipsis)))
}

prior.param <- function(prior.family, data.dim) {
  error.msg <- "Invalid argument to prior.family."

  if (!is.list(prior.family))
    stop(error.msg)

  # If a named list is found (rather than a list of named lists), it is
  #   interpreted as specifying the prior family (or families) for all factors.
  if (!is.null(names(prior.family[[1]])))
    prior.family <- list(prior.family)

  # Each top-level list element in prior.family corresponds to a factor.
  prior.family <- lapply(prior.family, function(k) {
    if (!is.list(k))
      stop(error.msg)
    if (length(k) == 1)
      k <- rep(k, data.dim)
    if (length(k) != data.dim)
      stop(error.msg)
    return(as.list(k))
  })

  prior.sign <- lapply(prior.family, lapply, `[[`, "sign")
  ebnm.fn    <- lapply(prior.family, lapply, `[[`, "ebnm.fn")

  return(list(prior.sign = prior.sign,
              ebnm.fn = ebnm.fn))
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
      # Default output columns for verbose.lvl = 3.
      param$verbose.fns       <- c(calc.obj.diff, calc.max.chg.EF)
      param$verbose.colnames  <- c("Obj Diff", "Max Chg")
      param$verbose.colwidths <- c(12, 12)
    } else if (verbose == -1) {
      # Default output columns for verbose.lvl = -1.
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

    if (length(chars) > 2) {
      stop("Unable to parse verbose output string.")
    } else if (length(chars) == 2) {
      n <- as.integer(chars[[2]])
      must.be.integer(n, lower = 1, upper = data.dim)
    } else {
      n <- NULL
    }

    if (chars[[1]] %in% c("O", "D")) {
      if (!is.null(n))
        warning("Dimension ignored for verbose objective output.")
      if (chars[[1]] == "O") {
        return(display.obj)
      } else { # if chars[[1]] == "D"
        return(calc.obj.diff)
      }
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

    if (chars[[1]] == "O") {
      return("Objective")
    } else if (chars[[1]] == "D") {
      return("Obj Diff")
    } else if (chars[[1]] == "L") {
      name <- "Max Chg"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    } else if (chars[[1]] == "W") {
      name <- "Whch"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    } else if (chars[[1]] == "S") {
      return(paste("Sparsity", chars[[2]]))
    } else if (chars[[1]] == "E") {
      return(paste("Excl", chars[[2]]))
    }

    stop("Unrecognized verbose output character.")
  })

  return(unlist(names))
}

look.up.verbose.colwidths <- function(verbose) {
  widths <- lapply(verbose, function(symbol) {
    char <- substring(symbol, 1, 1)

    if (char == "O") {
      return(16)
    } else if (char %in% c("D", "L", "S")) {
      return(13)
    } else if (char %in% c("W", "E")) {
      return(9)
    }

    stop("Unrecognized verbose output character.")
  })

  return(unlist(widths))
}
