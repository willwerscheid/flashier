# Sparse tensors. Internally, a list of sparse matrices (of class "Matrix"). It
# should be fastest to split up the matrices according to their smallest
# dimension (ie, for an m x n x p tensor with p < m and p < n, represent the
# tensor as a p-list of m x n matrices).
#
as.fullrank.sptensor <-function(x, ...) {
  return(sapply(x, as.matrix, simplify = "array"))
}

get.data.dims.sptensor <- function(x, ...) {
  return(c(dim(x[[1]]), length(x)))
}

get.data.dimnames.sptensor <- function(x, ...) {
  return(c(dimnames(x[[1]]), list(names(x))))
}

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

sq.nmode.prod.r1.sptensor <- function(x, r1, n, ...) {
  x2 <- lapply(x, `^`, 2)
  class(x2) <- "sptensor"
  return(nmode.prod.r1(x2, r1, n))
}
