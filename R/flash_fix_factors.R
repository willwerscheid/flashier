#' Fix flash factors
#'
#' Fixes or unfixes loadings \eqn{\ell_k} or factors \eqn{f_k} for one or more
#'   factor/loadings pairs. For a given pair, either the loadings or factor can
#'   be fixed --- but not both ---, and either all loadings or a subset can be
#'   fixed. To unfix, set \code{is.fixed = FALSE}.
#'
#' @param flash A \code{flash} or \code{flash.fit} object
#'
#' @param kset A vector of integers indexing the factor/loadings pairs whose
#'   loadings or factors are to be fixed (or unfixed).
#'
#' @param mode Set \code{mode = 1} to fix loadings and \code{mode = 2} to fix
#'   factors.
#'
#' @param is.fixed If \code{is.fixed = TRUE}, then all entries along the
#'   specified mode will be fixed. If only a subset are to be fixed,
#'   then \code{is.fixed} should be an appropriately-sized vector or
#'   matrix of values that can be coerced to logical. For example, if
#'   loadings for two factor/loadings pairs are to be fixed, then
#'   \code{is.fixed} can be a length-\eqn{n} vector or an \eqn{n} by 2
#'   matrix (where \eqn{n} is the number of rows in the data matrix).
#'   Finally, loadings can be "unfixed" by setting \code{is.fixed = FALSE}.
#'
#' @param use.fixed.in.ebnm By default, fixed elements are ignored when
#'   solving the EBNM subproblem in order to estimate the prior \eqn{\hat{g}}.
#'   This behavior can be changed by setting \code{use.fixed.in.ebnm = TRUE}.
#'   This is a global setting which applies to all factor/loadings pairs;
#'   behavior cannot vary from one pair to another.
#'
#' @export
#'
flash.fix.factors <- function(flash,
                              kset,
                              mode,
                              is.fixed = TRUE,
                              use.fixed.in.ebnm = NULL) {
  fit <- get.fit(flash)
  if (!is.null(use.fixed.in.ebnm)) {
    if (!is.logical(use.fixed.in.ebnm)) {
      stop("Argument to parameter use.fixed.in.ebnm must be logical.")
    }
    fit <- set.fixed.to.est.g(fit, use.fixed.in.ebnm)
  }

  must.be.valid.kset(fit, kset)
  must.be.integer(mode, lower = 1, upper = get.dim(fit), allow.null = FALSE)

  expect.nrow <- get.dims(fit)[mode]
  expect.ncol <- length(kset)

  if (!((is.vector(is.fixed) && length(is.fixed) %in% c(1, expect.nrow))
        || identical(dim(is.fixed), c(expect.nrow, expect.ncol)))) {
      stop("is.fixed must be a vector of length ", expect.nrow, " or a ",
           expect.nrow, " by ", expect.ncol, " matrix.")
  }
  is.fixed <- array(as.logical(is.fixed), dim = c(expect.nrow, expect.ncol))

  fix.dim <- get.fix.dim(fit)
  fix.idx <- get.fix.idx(fit)

  for (i in 1:length(kset)) {
    k <- kset[i]

    next.idx <- which(is.fixed[, i])

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
