#' Initialize a flash factor
#'
#' The default method for initializing the loadings \eqn{\ell_{\cdot k}} and
#'   factor values \eqn{f_{\cdot k}} of a new ("greedy") flash factor. It is
#'   essentially an implementation of the power method, but unlike many existing
#'   implementations, it can handle missing data and sign constraints. For details,
#'   see Chapter 2.2.3 in the reference below.
#'
#' @param flash A \code{flash_fit} object.
#'
#' @param sign_constraints This parameter can be used to constrain the sign of
#'   the initial factor and loadings. It should be a vector of length two with
#'   entries equal to -1, 0, or 1. The first entry constrains the sign of the
#'   loadings \eqn{\ell_{\cdot k}}, with -1 yielding nonpositive loadings, +1
#'   yielding nonnegative loadings, and 0 indicating that loadings should not be
#'   constrained. The second entry of \code{sign_constraints} similarly
#'   constrains the sign of factor values \eqn{f_{\cdot k}}. If
#'   \code{sign_constraints = NULL}, then no constraints will be applied.
#'
#' @param tol Convergence tolerance parameter. When the maximum (absolute)
#'   change over all values \eqn{\ell_{ik}} and \eqn{f_{jk}} is less than or
#'   equal to \code{tol}, initialization terminates. At each iteration, the
#'   factor and loadings are \eqn{L^2}-normalized. The default tolerance
#'   parameter is \eqn{\min(1 / n, 1 / p)}, where \eqn{n} is
#'   the number of rows in the data matrix and \eqn{p} is the number of columns.
#'
#' @param maxiter Maximum number of power iterations.
#'
#' @param seed Since initialization is random, a default seed is set for
#'   reproducibility.
#'
#' @return A list of length two consisting of, respectively, the vector of
#'   initial values for loadings \eqn{\ell_{\cdot k}} and the vector of initial
#'   factor values \eqn{f_{\cdot k}}.
#'
#' @seealso \code{\link{flash_greedy}},
#'   \code{\link{flash_greedy_init_softImpute}},
#'   \code{\link{flash_greedy_init_irlba}}
#'
#' @references
#' Jason Willwerscheid (2021), \emph{Empirical Bayes Matrix Factorization:
#'   Methods and Applications}. Ph.D. thesis, University of Chicago.
#'
#' @export
#'
flash_greedy_init_default <- function(flash,
                                      sign_constraints = NULL,
                                      tol = NULL,
                                      maxiter = 100,
                                      seed = 666) {
  set.seed(seed)
  if (is.null(sign_constraints)) {
    sign_constraints <- rep(0, get.dim(flash))
  }
  if (is.null(tol)) {
    tol <- 1 / max(get.dims(flash))
  }
  EF <- r1.random(get.dims(flash), sign_constraints)

  update.order <- 1:get.dim(flash)
  # Nonnegative/nonpositive dimensions are updated last.
  signed.dims <- which(sign_constraints %in% c(-1, 1))
  if (length(signed.dims) > 0) {
    which.signed <- which(update.order %in% signed.dims)
    update.order <- c(update.order[-which.signed], update.order[which.signed])
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.EF <- EF
    EF <- update.init.EF(EF, flash, update.order, sign_constraints)
    max.chg <- calc.max.abs.chg(EF, old.EF)
  }

  # Scale EF so that values aren't too different one dimension from another.
  EF <- scale.EF(EF)

  return(EF)
}

#' @exportS3Method NULL
update.init.EF <- function(EF, flash, update.order, sign_constraints) {
  if (is.null(sign_constraints))
    sign_constraints <- rep(0, get.dim(flash))

  for (n in update.order) {
    sign <- sign_constraints[n]
    EF <- update.init.EF.one.n(EF, n, flash, sign)
  }

  return(EF)
}

#' @exportS3Method NULL
update.init.EF.one.n <- function(EF, n, flash, sign) {
  R        <- get.R(flash)
  Y        <- get.Y(flash)
  Z        <- get.nonmissing(flash)
  flash.EF <- get.EF(flash)

  if (uses.R(flash)) {
    new.vals <- (nmode.prod.r1(R, EF[-n], n)
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  } else {
    new.vals <- ((nmode.prod.r1(Y, EF[-n], n)
                  - premult.nmode.prod.r1(Z, flash.EF, EF[-n], n))
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  }

  new.vals[is.na(new.vals)] <- 0

  if (sign == 1)
    new.vals <- pmax(new.vals, 0)
  if (sign == -1)
    new.vals <- pmin(new.vals, 0)

  EF[[n]] <- new.vals

  return(EF)
}

#' @exportS3Method NULL
scale.EF <- function(EF) {
  norms <- lapply(EF, function(x) {sqrt(sum(x^2))})

  if (all(unlist(norms) > 0)) {
    EF <- mapply(`/`, EF, norms, SIMPLIFY = FALSE)
    EF <- lapply(EF, `*`, prod(unlist(norms))^(1/length(EF)))
  } else {
    warning("Fitting stopped after the initialization function failed to find",
            " a non-zero factor.")
    EF <- lapply(EF, `*`, 0)
  }
  class(EF) <- "r1"

  return(EF)
}

#' Initialize a flash factor using softImpute
#'
#' Initializes a new ("greedy") flash factor using \code{\link[softImpute]{softImpute}}.
#'   When there is missing data, this can yield better results than
#'   \code{\link{flash_greedy_init_default}} without sacrificing much (if any) speed.
#'
#' @inheritParams flash_greedy_init_default
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link[softImpute]{softImpute}}.
#'
#' @inherit flash_greedy_init_default return
#'
#' @seealso \code{\link{flash_greedy}},
#'   \code{\link{flash_greedy_init_default}},
#'   \code{\link{flash_greedy_init_irlba}}
#'
#' @importFrom stats residuals
#' @importFrom softImpute softImpute
#'
#' @export
#'
flash_greedy_init_softImpute <- function(flash, seed = 666, ...) {
  set.seed(seed)

  if (get.dim(flash) > 2)
    stop("softImpute cannot be used with tensors.")

  if (inherits(get.Y(flash), "lowrank"))
    stop("softImpute cannot be used with low-rank matrix representations.")

  if (inherits(get.Y(flash), "lrps"))
    stop("softImpute cannot be used with low-rank plus sparse representations.")

  si.res <- softImpute(residuals(flash), rank.max = 1, ...)
  EF <- list(si.res$u * sqrt(si.res$d), si.res$v * sqrt(si.res$d))

  return(EF)
}

#' Initialize a flash factor using IRLBA
#'
#' Initializes a new ("greedy") flash factor using \code{\link[irlba]{irlba}}. This
#'   can be somewhat faster than \code{\link{flash_greedy_init_default}} for large,
#'   dense data matrices. For sparse matrices of class \code{Matrix}, the
#'   default initialization should generally be preferred.
#'
#' @inheritParams flash_greedy_init_default
#'
#' @param ... Additional parameters to be passed to \code{\link[irlba]{irlba}}.
#'
#' @inherit flash_greedy_init_default return
#'
#' @seealso \code{\link{flash_greedy}},
#'   \code{\link{flash_greedy_init_default}},
#'   \code{\link{flash_greedy_init_softImpute}}
#'
#' @importFrom stats residuals
#' @importFrom irlba irlba
#'
#' @export
#'
flash_greedy_init_irlba <- function(flash, seed = 666, ...) {
  set.seed(seed)

  if (get.dim(flash) > 2)
    stop("irlba cannot be used with tensors.")

  if (inherits(get.Y(flash), "lowrank"))
    stop("irlba cannot be used with low-rank matrix representations.")

  if (any_missing(flash)) {
    stop("irlba cannot be used when there is missing data.")
  }

  irlba.res <- irlba(residuals(flash), nv = 1, nu = 1, ...)

  EF <- list(
    as.vector(irlba.res$u * sqrt(irlba.res$d)),
    as.vector(irlba.res$v * sqrt(irlba.res$d))
  )

  return(EF)
}
