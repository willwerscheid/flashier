#' Fitted method for flash objects
#'
#' Given a \code{\link{flash}} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @return The matrix of "fitted values."
#'
#' @importFrom stats fitted
#' @method fitted flash
#'
#' @export
#'
fitted.flash <- function(object, ...) {
  return(fitted.flash_fit(get.fit(object)))
}

#' Fitted method for flash fit objects
#'
#' Given a \code{\link{flash_fit}} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash_fit}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @return The matrix of "fitted values."
#'
#' @importFrom stats fitted
#' @method fitted flash_fit
#'
#' @export
#'
fitted.flash_fit <- function(object, ...) {
  if (get.n.factors(object) == 0) {
    warning("Flash object does not have any factors.")
    return(matrix(data = 0, nrow = nrow(object$Y), ncol = ncol(object$Y)))
  }
  if (get.dim(object) > 2) {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }

  return(do.call(tcrossprod, get.EF(object)))
}

#' Residuals method for flash objects
#'
#' Given a \code{\link{flash}} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash
#'
#' @return The matrix of expected residuals.
#'
#' @importFrom stats residuals
#' @method residuals flash
#'
#' @export
#'
residuals.flash <- function(object, ...) {
  return(residuals.flash_fit(get.fit(object)))
}

#' Residuals method for flash fit objects
#'
#' Given a \code{\link{flash_fit}} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash_fit
#'
#' @return The matrix of expected residuals.
#'
#' @importFrom stats residuals
#' @method residuals flash_fit
#'
#' @export
#'
residuals.flash_fit <- function(object, ...) {
  if (uses.R(object)) {
    R <- get.R(object)
  } else {
    R <- get.Y(object) - lowrank.expand(get.EF(object))
  }

  if (any_missing(object)) {
    R[get.nonmissing(object) == 0] <- NA
  }

  return(R)
}

#' LDF method for flash and flash fit objects
#'
#' Given a \code{\link{flash}} or \code{\link{flash_fit}} object, returns the LDF
#'   decomposition \eqn{Y \approx LDF'}.
#'
#' When the prior families \eqn{G_\ell^{(k)}} and \eqn{G_f^{(k)}} are closed
#'   under scaling (as is typically the case), then the EBMF model (as
#'   described in the documention to function \code{\link{flash}}) is only
#'   identifiable up to scaling by a diagonal matrix \eqn{D}:
#'   \deqn{Y = LDF' + E.}
#'
#' Method \code{ldf} scales columns \eqn{\ell_k} and \eqn{f_k}
#'   so that, depending on the argument to parameter \code{type}, their
#'   1-norms, 2-norms, or infinity norms are equal to 1.
#'
#' @param object An object inheriting from class \code{flash} or
#'   \code{flash_fit}.
#'
#' @param type Takes identical arguments to function \code{\link[base]{norm}}. Use
#'   \code{"f"} or \code{"2"} for the 2-norm (Euclidean norm); \code{"o"} or
#'   \code{"1"} for the 1-norm (taxicab norm); and \code{"i"} or \code{"m"} for
#'   the infinity norm (maximum norm).
#'
#' @return A list with fields \code{L}, \code{D}, and \code{F}, each of which
#'   corresponds to one of the matrices in the decomposition \eqn{Y \approx LDF'}
#'   (with the columns of \eqn{L} and \eqn{F} scaled according to
#'   argument \code{type}). Note that \code{D} is returned as a vector rather
#'   than a matrix (the vector of diagonal entries in \eqn{D}). Thus, "fitted
#'   values" \eqn{LDF'} can be recovered as \code{L \%*\% diag(D) \%*\% t(F)}.
#'
#' @export
#'
ldf <- function(object, type) {
  UseMethod("ldf", object)
}

#' @describeIn ldf LDF decomposition for \code{\link{flash}} objects
#'
#' @method ldf flash
#'
#' @export
#'
ldf.flash <- function(object, type = "f") {
  return(ldf.flash_fit(get.fit(object), type = type))
}

#' @describeIn ldf LDF decomposition for \code{\link{flash_fit}} objects
#'
#' @method ldf flash_fit
#'
#' @export
#'
ldf.flash_fit <- function(object, type = "f") {
  type <- tolower(type)
  if (type == "2") {
    type <- "f"
  } else if (type == "1") {
    type <- "o"
  } else if (type == "m") {
    type <- "i"
  }

  if (get.n.factors(object) == 0) {
    stop("Flash fit does not have any factors.")
  }
  if (get.dim(object) > 2) {
    stop("S3 method \"ldf\" not available for tensors.")
  }

  ldf <- calc.normalized.loadings(object, type = type)

  ret <- list()
  ret$L <- ldf$normalized.loadings[[1]]
  ret$D <- ldf$scale.constants
  ret$F <- ldf$normalized.loadings[[2]]

  return(ret)
}

calc.normalized.loadings <- function(flash, for.pve = FALSE, type = "f") {
  ret <- list()

  if (for.pve) {
    loadings <- get.EF2(flash)
    norms <- lapply(loadings, colSums)
  } else {
    loadings <- get.EF(flash)
    if (type == "f") {
      norms <- lapply(loadings, function(x) {sqrt(colSums(x^2))})
    } else if (type == "o") {
      norms <- lapply(loadings, function(x) {colSums(abs(x))})
    } else if (type == "i") {
      norms <- lapply(loadings, function(x) {apply(abs(x), 2, max)})
    } else {
      stop("Norm type not recognized.")
    }
  }

  # Zero factors are "normalized" to zero.
  norms <- lapply(norms, function(x) {x[is.zero(flash)] <- Inf; x})
  L <- mapply(loadings, norms, FUN = function(X, y) {
    X / rep(y, each = nrow(X))
  }, SIMPLIFY = FALSE)
  # if (!for.pve) {
  #   L2 <- mapply(get.EF2(flash), norms, FUN = function(X, y) {
  #     X / rep(y^2, each = nrow(X))
  #   }, SIMPLIFY = FALSE)
  #   SD <- mapply(L2, L,
  #                FUN = function(EX2, EX) {sqrt(pmax(EX2 - EX^2, 0))},
  #                SIMPLIFY = FALSE)
  # }

  L <- propagate.names(L, flash)

  norms <- do.call(rbind, norms)
  ret$scale.constants <- apply(norms, 2, prod)
  ret$scale.constants[is.zero(flash)] <- 0
  ret$normalized.loadings <- L
  # if (!for.pve)
  #   ret$loading.SDs <- SD

  return(ret)
}

propagate.names <- function(mats, flash) {
  data.dimnames <- get.dimnames(flash)
  for (n in 1:get.dim(flash)) {
    if (!is.null(mats[[n]]) && !is.null(data.dimnames) && !is.null(data.dimnames[[n]]))
      rownames(mats[[n]]) <- data.dimnames[[n]]
  }
  return(mats)
}

#' @export
print.flash = function(x, ...) {
  if (x$n_factors == 0) {
    cat("Flash object with zero factors.\n")
  } else if (x$n_factors == 1) {
    cat("Flash object with one factor.\n")
    cat(sprintf("  Proportion of variance explained: %0.3f\n", x$pve))
  } else {
    cat(sprintf("Flash object with %d factors.\n", x$n_factors))
    if (all(x$pve < 0.001)) {
      cat("  All factors have PVE < 0.001.")
    } else {
      cat("  Proportion of variance explained")
      if (any(x$pve < 0.001)) {
        cat("*")
      }
      cat(":\n")
      for (k in 1:length(x$pve)) {
        if (x$pve[k] > 0.001) {
          cat(sprintf("    Factor %d: %0.3f\n", k, x$pve[k]))
        }
      }
      if (any(x$pve < 0.001)) {
        cat(paste("    *Factors with PVE < 0.001 are omitted from this",
                  "summary.\n"))
      }
    }
  }

  if (!is.na(x$elbo)) {
    cat(sprintf("  Variational lower bound: %0.3f\n", x$elbo))
  }
}
