#' Initialize a flash factor
#'
#' The default method for initializing a new flash factor.
#'
#' @param flash A \code{flash.fit} object.
#'
#' @param dim.signs This parameter can be used to constrain the sign of the
#'   initial loadings. It should be a vector of length two with entries equal
#'   to -1, 0, or 1. The first entry dictates the sign of the loadings
#'   \eqn{\ell_k}, with -1 yielding nonpositive loadings, +1 yielding
#'   nonnegative loadings, and 0 indicating that loadings should not be
#'   constrained. The second entry of \code{dim.signs} similarly constrains
#'   the sign of the factor \eqn{f_k}.
#'
#' @param tol Convergence tolerance.
#'
#' @param maxiter Maximum number of iterations.
#'
#' @param seed Since initialization is random, a default seed is set for
#'   reproducibility.
#'
#' @seealso \code{\link{init.fn.softImpute}}, \code{\link{init.fn.irlba}}
#'
#' @export
#'
init.fn.default <- function(flash,
                            dim.signs = rep(0, get.dim(flash)),
                            tol = 1 / max(get.dims(flash)),
                            maxiter = 100,
                            seed = 666) {
  set.seed(seed)
  EF <- r1.random(get.dims(flash), dim.signs)

  update.order <- 1:get.dim(flash)
  # Nonnegative/nonpositive dimensions are updated last.
  signed.dims <- which(dim.signs %in% c(-1, 1))
  if (length(signed.dims) > 0) {
    which.signed <- which(update.order %in% signed.dims)
    update.order <- c(update.order[-which.signed], update.order[which.signed])
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.EF <- EF
    EF <- update.init.EF(EF, flash, update.order, dim.signs)
    max.chg <- calc.max.abs.chg(EF, old.EF)
  }

  # Scale EF so that values aren't too different one dimension from another.
  EF <- scale.EF(EF)

  return(EF)
}

update.init.EF <- function(EF, flash, update.order, dim.signs) {
  if (is.null(dim.signs))
    dim.signs <- rep(0, get.dim(flash))

  for (n in update.order) {
    sign <- dim.signs[n]
    EF <- update.init.EF.one.n(EF, n, flash, sign)
  }

  return(EF)
}

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
#' Initializes a new flash factor using \code{\link[softImpute]{softImpute}}.
#'   When there is missing data, this can yield better results than
#'   \code{\link{init.fn.default}} without sacrificing much (if any) speed.
#'
#' @inheritParams init.fn.default
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link[softImpute]{softImpute}}.
#'
#' @seealso \code{\link{init.fn.default}}, \code{\link{init.fn.irlba}}
#'
#' @importFrom softImpute softImpute
#'
#' @export
#'
init.fn.softImpute <- function(flash, seed = 666, ...) {
  set.seed(seed)

  if (get.dim(flash) > 2)
    stop("softImpute cannot be used with tensors.")

  if (inherits(get.Y(flash), "lowrank"))
    stop("softImpute cannot be used with low-rank matrix representations.")

  si.res <- softImpute(fitted(flash), rank.max = 1, ...)
  EF <- list(si.res$u * sqrt(si.res$d), si.res$v * sqrt(si.res$d))

  return(EF)
}

#' Initialize a flash factor using IRLBA
#'
#' Initializes a new flash factor using \code{\link[irlba]{irlba}}. This
#'   can be somewhat faster than \code{\link{init.fn.default}} for large,
#'   dense data matrices. For sparse matrices of class \code{Matrix}, the
#'   default initialization should generally be preferred.
#'
#' @inheritParams init.fn.default
#'
#' @param ... Additional parameters to be passed to \code{\link[irlba]{irlba}}.
#'
#' @seealso \code{\link{init.fn.default}}, \code{\link{init.fn.softImpute}}
#'
#' @importFrom irlba irlba
#'
#' @export
#'
init.fn.irlba <- function(flash, seed = 666, ...) {
  set.seed(seed)

  if (get.dim(flash) > 2)
    stop("irlba cannot be used with tensors.")

  if (inherits(get.Y(flash), "lowrank"))
    stop("irlba cannot be used with low-rank matrix representations.")

  if (any.missing(flash)) {
    stop("irlba cannot be used when there is missing data.")
  }

  irlba.res <- irlba(fitted(flash), nv = 1, nu = 1, ...)

  EF <- list(
    as.vector(irlba.res$u * sqrt(irlba.res$d)),
    as.vector(irlba.res$v * sqrt(irlba.res$d))
  )

  return(EF)
}
