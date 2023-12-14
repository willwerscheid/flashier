# Sparse tensors. Internally, a list of sparse matrices (of class "Matrix"). It
# should be fastest to split up the matrices according to their smallest
# dimension (ie, for an m x n x p tensor with p < m and p < n, represent the
# tensor as a p-list of m x n matrices).
#
#' @exportS3Method as.fullrank sptensor
as.fullrank.sptensor <-function(x, ...) {
  return(sapply(x, as.matrix, simplify = "array"))
}

#' @exportS3Method get.data.dims sptensor
get.data.dims.sptensor <- function(x, ...) {
  return(c(dim(x[[1]]), length(x)))
}

#' @exportS3Method get.data.dimnames sptensor
get.data.dimnames.sptensor <- function(x, ...) {
  return(c(dimnames(x[[1]]), list(names(x))))
}

#' @exportS3Method nmode.prod.vec sptensor
nmode.prod.vec.sptensor <- function(x, v, n, ...) {
  if (n == 1 && identical(v, 1))
    return(Reduce(cbind, lapply(x, colSums)))
  if (n == 2 && identical(v, 1))
    return(Reduce(cbind, lapply(x, rowSums)))
  if (n == 1 || n == 2)
    return(Reduce(cbind, lapply(x, nmode.prod.vec, v = v, n = n)))
  if (n == 3 && identical(v, 1))
    return(Reduce(`+`, x))
  if (n == 3) {
    return(Reduce(`+`, mapply(`*`, x, v)))
  }
}

#' @exportS3Method sq.nmode.prod.r1 sptensor
sq.nmode.prod.r1.sptensor <- function(x, r1, n, ...) {
  x2 <- lapply(x, `^`, 2)
  class(x2) <- "sptensor"
  return(nmode.prod.r1(x2, r1, n))
}
