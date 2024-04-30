#   N-mode products generalize matrix-vector products and constitute the key
# operations for this package. They currently work for matrices, sparse
# matrices (as implemented by package Matrix), 3-dimensional tensors (which
# should be input as arrays), low-rank matrices and tensors (objects of
# class "lowrank" which are stored as lists of matrices of dimension d_n x k),
# and "low-rank plus sparse" matrix representations (see data_lrps.R).
# Ideally, they would also accept sparse tensors and tensors of arbitrary
# dimension, but I am not aware of any good R packages for handling these
# types of object.
#
#   To implement new data structures, each of the following methods are
# required:
#
# REQUIRED -----

as.fullrank <- function(x, ...) {
  UseMethod("as.fullrank", x)
}

get.data.dims <- function(x, ...) {
  UseMethod("get.data.dims", x)
}

get.data.dimnames <- function(x, ...) {
  UseMethod("get.data.dimnames", x)
}

# Takes three arguments: X, v, and n.
# X is a m_1 x m_2 (x m_3) matrix (array) and v is an m_n-vector. The n-mode
#   product multiplies the jth n-slice of X by v_j and then sums the slices
#   together. For example, the 3-mode product of a 3-dimensional tensor is
#   \sum_k v_k X_{ijk}. v can be set to 1 to indicate an m_n-vector of all
#   ones.
#
nmode.prod.vec <- function(x, ...) {
  UseMethod("nmode.prod.vec", x)
}

# Takes three arguments: X, r1, and n.
# X is a m_1 x m_2 (x m_3) matrix (array), r1 is a r1 object whose dimension
#   is one less than X, and n indicates the dimension that is omitted. For
#   matrices, this is simply the usual matrix product of X and the vector
#   in r1, so the default method will usually suffice for two-dimensional
#   data structures. For tensors, this is the n-mode product of X and the
#   rank-one matrix formed by the outer product of the vectors in r1. For
#   example, if r1 = list(u, v) and n = 1, then \sum_{j, k} u_j v_k X_{ijk}
#   is returned. Any of the list elements in r1 can be set to 1 to indicate
#   a vector of all ones.
#
nmode.prod.r1 <- function(x, ...) {
  UseMethod("nmode.prod.r1", x)
}

# Takes three arguments: X, r1, and n.
# Calculates the n-mode product between X^2 and an r1 object.
sq.nmode.prod.r1 <- function(x, ...) {
  UseMethod("sq.nmode.prod.r1", x)
}

# Takes four arguments: X, lr, r1, and n.
# The following calculates the n-mode product between the matrix or tensor
#   X * LF (where LF is stored as the "lowrank" object lr and multiplication
#   is elementwise) and an r1 object. The idea is both to avoid forming a large
#   matrix (or tensor) and to be able to take advantage of any sparsity in X.
#   Methods should also be able to handle the case where lr is a full-rank
#   matrix or tensor. (TODO: Eliminate need to handle full-rank objects.)
#
premult.nmode.prod.r1 <- function(x, ...) {
  UseMethod("premult.nmode.prod.r1", x)
}

# NOT USED -----

# Only needed for parallel backfits (which are not currently available):
#
nmode.prod.lowrank <- function(x, ...) {
  UseMethod("nmode.prod.lowrank", x)
}

# DEFAULT METHODS -----

#' @export
as.fullrank.default <- function(x, ...) {
  return(x)
}

#' @export
get.data.dims.default <- function(x, ...) {
  return(dim(x))
}

#' @export
get.data.dimnames.default <- function(x, ...) {
  return(dimnames(x))
}

nmode.ops.error <- paste("N-mode products are not yet implemented for",
                         "that data structure.")

# Imports needed for sparse matrices:
#' @importFrom Matrix colSums rowSums t crossprod tcrossprod
#' @export
nmode.prod.vec.default <- function(x, v, n, ...) {
  if (length(get.data.dims(x)) == 2) {
    if (n == 1 && identical(v, 1))
      return(colSums(x))
    if (n == 1)
      return(as.vector(v %*% x))
    if (n == 2 && identical(v, 1))
      return(rowSums(x))
    if (n == 2)
      return(as.vector(x %*% v))
  } else if (length(get.data.dims(x)) == 3) {
    if (n == 1 && identical(v, 1))
      return(apply(x, 3, colSums))
    if (n == 1)
      return(apply(x, 3, FUN = function(M) {t(v) %*% M}))
    if (n == 2 && identical(v, 1))
      return(apply(x, 3, rowSums))
    if (n == 2)
      return(apply(x, 3, FUN = function(M) {M %*% v}))
    if (n == 3 && identical(v, 1))
      return(apply(x, 2, rowSums))
    if (n == 3)
      return(apply(x, 2, FUN = function(M) {M %*% v}))
  } else {
    stop(nmode.ops.error)
  }
}

#' @export
nmode.prod.r1.default <- function(x, r1, n, ...) {
  if (is.null(x)) {
    return(0)
  } else if (identical(x, 1)) {
    return(r1.sum(r1))
  } else if (length(get.data.dims(x)) == 2) {
    return(nmode.prod.vec(x, unlist(r1), (1:2)[-n]))
  } else if (length(get.data.dims(x)) == 3) {
    ns <- (1:3)[-n]
    # Go backwards (otherwise indices get messed up).
    return(nmode.prod.vec(nmode.prod.vec(x, r1[[2]], ns[2]), r1[[1]], ns[1]))
  } else {
    stop(nmode.ops.error)
  }
}

#' @export
sq.nmode.prod.r1.default <- function(x, r1, n, ...) {
  return(nmode.prod.r1(x^2, r1, n))
}

#' @export
premult.nmode.prod.r1.default <- function(x, lr, r1, n, ...) {
  if (is.null(lr))
    return(0)

  if (identical(x, 1)) {
    return(nmode.prod.r1(lr, r1, n))
  }

  if (inherits(lr, "lowrank")) {
    if (length(lr) == 2) {
      u <- unlist(r1) * lr[[-n]]
      if (n == 1)
        return(rowSums(lr[[1]] * (x %*% u)))
      if (n == 2)
        return(colSums(t(lr[[2]]) * (t(u) %*% x)))
    } else if (length(lr) == 3) {
      ns <- (1:3)[-n]
      # Sum over third (or second) index. Get k m_1 x m_2 (or m_1 x m_3) matrices
      #   stacked on top of one another.
      tmp <- nmode.prod.vec(x, r1[[2]] * lr[[ns[2]]], ns[[2]])
      if (inherits(tmp, "lowrank"))
        tmp <- lowrank.expand(tmp)
      # Elementwise multiply by an m_1 x k or m_2 x k matrix, then sum over the
      #   m_ index.
      if (ns[[1]] == 1) {
        tmp <- tmp * as.vector(lr[[1]] * r1[[1]])
        tmp <- colSums(matrix(tmp, nrow = dim(x)[1]))
        # Have m_n k vectors concatenated one after another.
        return(rowSums(lr[[n]]
                       * matrix(tmp, nrow = dim(x)[n], byrow = TRUE)))
      } else { # ns[[1]] == 2
        tmp <- tmp * rep(t(lr[[2]] * r1[[1]]), each = dim(x)[1])
        tmp <- rowSums(tmp)
        # Have k m_n vectors concatenated one after another.
        return(rowSums(lr[[n]] * matrix(tmp, nrow = dim(x)[n])))
      }
    } else {
      stop(nmode.ops.error)
    }
  }

  # else lr is full-rank:
  return(nmode.prod.r1(x * lr, r1, n))
}

# Not currently used -----

# Only needed for parallel backfits:
#' @export
nmode.prod.lowrank.default <- function(x, Y, n, ...) {
  if (length(get.data.dims(x)) == 2) {
    if (n == 1)
      return(as.matrix(x %*% Y[[1]]))
    if (n == 2)
      return(t(as.matrix(t(Y[[1]]) %*% x)))
  } else if (length(get.data.dims(x)) == 3) {
    stop("N-mode products of tensors with matrices have not yet been ",
         "implemented.")
  } else {
    stop(nmode.ops.error)
  }
}
