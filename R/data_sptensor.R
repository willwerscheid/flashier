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
  # cbind names the first column, which can throw weird errors.
  cbind0 <- function(...) cbind(..., deparse.level = 0)
  if (n == 1 && identical(v, 1))
    return(Reduce(cbind0, lapply(x, colSums)))
  if (n == 2 && identical(v, 1))
    return(Reduce(cbind0, lapply(x, rowSums)))
  if (n == 1 || n == 2)
    return(Reduce(cbind0, lapply(x, nmode.prod.vec, v = v, n = n)))
  if (n == 3 && identical(v, 1))
    return(Reduce(`+`, x))
  if (n == 3) {
    return(Reduce(`+`, mapply(`*`, x, v)))
  }
}

#' @exportS3Method premult.nmode.prod.r1 sptensor
premult.nmode.prod.r1.sptensor <- function(x, lr, r1, n, ...) {
  if (is.null(lr))
    return(0)

  if (inherits(lr, "lowrank")) {
    ns <- (1:3)[-n]
    n1 <- ns[1]
    n2 <- ns[2]
    K <- ncol(lr[[1]])

    # Sum over third (or second) index. Get a length-K list of m_1 x m_2
    #   (or m_1 x m_3) sparse matrices.
    tmp <- apply(lr[[n2]], 2, function(lr_k) {
      nmode.prod.vec(x, r1[[2]] * lr_k, n2)
    }, simplify = FALSE)
    # Multiply each matrix by an m_1 (or m_2) vector using broadcasting and
    #   take sums; gives m_2/m_3 x K (or m_1 x K) matrix.
    if (n1 == 1) {
      tmp <- sapply(1:K, function(k) {
        colSums(tmp[[k]] * lr[[n1]][, k] * r1[[1]])
      })
    } else { # n1 == 2, so n2 == 3
      tmp <- sapply(1:K, function(k) {
        colSums(t(tmp[[k]]) * lr[[n1]][, k] * r1[[1]])
      })
    }
    return(rowSums(tmp * lr[[n]]))
  }

  # else lr is full-rank (TODO: shouldn't get here):
  return(nmode.prod.r1(as.fullrank(x) * lr, r1, n))
}

#' @exportS3Method sq.nmode.prod.r1 sptensor
sq.nmode.prod.r1.sptensor <- function(x, r1, n, ...) {
  x2 <- lapply(x, `^`, 2)
  class(x2) <- "sptensor"
  return(nmode.prod.r1(x2, r1, n))
}
