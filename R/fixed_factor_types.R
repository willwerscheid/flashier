#' Flashier fixed factors
#'
#' Functions for creating different types of fixed factors. Note that loadings
#'   can only be fixed along one mode (e.g., rows, columns) at a time.
#'   \code{ones.factor} creates a factor whose loadings along one mode are all
#'   fixed at 1. This type of factor can be used (for example) to estimate
#'   row and column means. \code{sparse.factors} creates one or more factors
#'   with some loadings fixed at zero (and with all other loadings to be
#'   estimated). \code{fixed.factors} creates one or more factors with loadings
#'   fixed at arbitrarily specified values.
#'
#' @param n The mode of the fixed loadings. For matrices, \code{n = 1}
#'   indicates that row loadings will be fixed. \code{n = 2} designates
#'   columns.
#'
#' @rdname fixed.factor
#'
#' @export
#'
ones.factor <- function(n) {
  return(list(list(dim = n, vals = 1)))
}

#' @param nz.idx For sparse factors, \code{nz.idx} gives the indices of the
#'   nonzero loadings (that is, the loadings that are to be estimated rather
#'   than fixed at zero). If \code{nz.idx} is a vector, then a single sparse
#'   factor will be created. If it is a list of vectors, then multiple factors
#'   will be created.
#'
#' @rdname fixed.factor
#'
#' @export
#'
sparse.factors <- function(n, nz.idx) {
  if (!is.list(nz.idx))
    nz.idx <- list(nz.idx)
  ff <- lapply(nz.idx, function(elem) {
    list(dim = n, vals = 0, idx = -elem)
  })
  return(ff)
}

#' @param vals The values of the fixed loadings. Can be a vector, which is
#'   interpreted as a single fixed factor, or a matrix, each column of which
#'   is interpreted as a different fixed factor.
#'
#' @param idx The indices of the fixed loadings. Set \code{idx = NULL} to fix
#'   all mode-\code{n} loadings.
#'
#' @rdname fixed.factor
#'
#' @export
#'
fixed.factors <- function(n, vals, idx = NULL) {
  if (is.matrix(vals)) {
    ff <- apply(vals, 2, function(col) {
      list(dim = n, vals = col, idx = idx)
    })
  } else {
    ff <- list(list(dim = n, vals = vals, idx = idx))
  }
  return(ff)
}
