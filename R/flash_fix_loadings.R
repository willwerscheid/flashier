#' Fix flash loadings
#'
#' Fixes some or all loadings within one or more flash factors.
#'
#' @param flash A \code{flash} or \code{flash.fit} object
#'
#' @param kset A vector of integers indexing the factors whose loadings are to
#'   be fixed (or unfixed).
#'
#' @param mode The mode along which loadings are to be fixed (for matrix data,
#'   \code{mode = 1} fixes row loadings and \code{mode = 2} fixes column
#'   loadings). For any given factor, loadings may only be fixed along a single
#'   mode.
#'
#' @param is.fixed If \code{is.fixed = TRUE}, then all loadings along the
#'   specified mode will be fixed. If only a subset of loadings are to be
#'   fixed, then \code{is.fixed} should be an appropriately-sized vector or
#'   matrix of values that can be coerced to logical. For example, if row
#'   loadings for two factors are to be fixed, then \code{is.fixed} can be
#'   a length-n vector or an n by 2 matrix (where n is the number of rows in
#'   the data matrix). Finally, loadings can be "unfixed" by setting
#'   \code{is.fixed = FALSE}.
#'
#' @export
#'
flash.fix.loadings <- function(flash, kset, mode, is.fixed = TRUE) {
  fit <- get.fit(flash)

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

    if (length(fix.dim) >= k && !is.null(fix.dim[[k]]) && fix.dim[[k]] != mode) {
        stop(paste("Loadings can only be fixed along a single mode for any",
                   "given factor."))
    }

    fix.dim[[k]] <- mode
    fix.idx[[k]] <- which(is.fixed[, i])
  }

  fit <- set.fix.dim(fit, fix.dim)
  fit <- set.fix.idx(fit, fix.idx)

  flash <- set.fit(flash, fit)

  return(flash)
}
