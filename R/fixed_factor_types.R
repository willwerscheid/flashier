#' @export
#'
mean.factor <- function(n) {
  return(list(list(dim = n, vals = 1)))
}

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
