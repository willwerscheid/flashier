#' Fix flash factors
#'
#' Fixes loadings \eqn{\ell_{\cdot k}} or factors \eqn{f_{\cdot k}}
#'   for one or more factor/loadings pairs, so that their values are not
#'   updated during subsequent backfits. For a given pair, either the loadings
#'   or factor can be fixed, but not both, and either all entries or a subset
#'   can be fixed. To unfix, use function \code{\link{flash_factors_unfix}}. See
#'   \code{\link{flash_factors_init}} for an example of usage.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param kset A vector of integers indexing the factor/loadings pairs whose
#'   loadings or factors are to be fixed.
#'
#' @param which_dim Whether to fix factors or loadings.
#'
#' @param fixed_idx If \code{fixed_idx = NULL}, then all loadings or factor
#'   values will be fixed. If only a subset are to be fixed,
#'   then \code{fixed_idx} should be an appropriately-sized vector or
#'   matrix of values that can be coerced to logical. For example, if
#'   a subset of loadings for two factor/loadings pairs are to be fixed, then
#'   \code{fixed_idx} should be a length-\eqn{n} vector or an \eqn{n} by 2
#'   matrix (where \eqn{n} is the number of rows in the data matrix).
#'
#' @param use_fixed_in_ebnm By default, fixed elements are ignored when
#'   solving the EBNM subproblem in order to estimate the prior \eqn{\hat{g}}.
#'   This behavior can be changed by setting \code{use_fixed_in_ebnm = TRUE}.
#'   This is a global setting which applies to all factor/loadings pairs;
#'   behavior cannot vary from one factor/loadings pair to another.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   factors or loadings fixed as specified.
#'
#' @export
#'
flash_factors_fix <- function(flash,
                              kset,
                              which_dim = c("factors", "loadings"),
                              fixed_idx = NULL,
                              use_fixed_in_ebnm = NULL) {
  fit <- get.fit(flash)
  if (!is.null(use_fixed_in_ebnm)) {
    if (!is.logical(use_fixed_in_ebnm)) {
      stop("Argument to parameter use_fixed_in_ebnm must be logical.")
    }
    fit <- set.fixed.to.est.g(fit, use_fixed_in_ebnm)
  }

  must.be.valid.kset(fit, kset)
  which_dim <- match.arg(which_dim)
  if (which_dim == "factors") {
    mode <- 2
  } else { # which_dim == "loadings"
    mode <- 1
  }

  expect.nrow <- get.dims(fit)[mode]
  expect.ncol <- length(kset)

  if (is.null(fixed_idx)) {
    fixed_idx <- TRUE
  }
  if (!((is.vector(fixed_idx) && length(fixed_idx) %in% c(1, expect.nrow))
        || identical(dim(fixed_idx), c(expect.nrow, expect.ncol)))) {
    stop("fixed_idx must be a vector of length ", expect.nrow, " or a ",
         expect.nrow, " by ", expect.ncol, " matrix.")
  }
  fixed_idx <- array(as.logical(fixed_idx), dim = c(expect.nrow, expect.ncol))

  fix.dim <- get.fix.dim(fit)
  fix.idx <- get.fix.idx(fit)

  for (i in 1:length(kset)) {
    k <- kset[i]

    next.idx <- which(fixed_idx[, i])

    if (length(next.idx) == 0) {
      fix.dim[[k]] <- mode
      fix.idx[[k]] <- numeric(0)
    } else if (length(fix.dim) >= k
               && !is.null(fix.dim[[k]])
               && fix.dim[[k]] != mode) {
      stop(paste("Loadings can only be fixed along a single mode for any",
                 "given factor."))
    } else {
      fix.dim[[k]] <- mode
      fix.idx[[k]] <- next.idx
    }
  }

  fit <- set.fix.dim(fit, fix.dim)
  fit <- set.fix.idx(fit, fix.idx)

  flash <- set.fit(flash, fit)

  return(flash)
}

#' Unfix flash factors
#'
#' If loadings \eqn{\ell_{\cdot k}} or factors \eqn{f_{\cdot k}} for one or
#'   more factor/loadings pairs have been "fixed" using function
#'   \code{\link{flash_factors_fix}}, then they can be unfixed using
#'   function \code{flash_factors_unfix}.
#'
#' @inheritParams flash_factors_fix
#'
#' @param kset A vector of integers indexing the factor/loadings pairs whose
#'   values are to be unfixed.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   values for the factor/loadings pairs specified by \code{kset} unfixed.
#'
#' @export
#'
flash_factors_unfix <- function(flash, kset) {
  return(flash_factors_fix(flash, kset, fixed_idx = FALSE))
}
