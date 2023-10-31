#   N-mode products generalize matrix-vector products and constitute the key
# operations for this package. They currently work for matrices, sparse
# matrices (as implemented by package Matrix), 3-dimensional tensors (which
# should be input as arrays), and low-rank matrices and tensors (objects of
# class "lowrank" which are stored as lists of matrices of dimension d_n x k).
# Ideally, they would also accept sparse tensors and tensors of arbitrary
# dimension, but I am not aware of any good R packages for handling these types
# of object.
#
# Added Nov 2023: low-rank plus sparse matrix representations (objects of
# class "lrps" which are stored as lists of length two with field "lowrank"
# containing an object of class "lowrank" and field "sparse" containing a
# matrix, ideally but not necessarily of class Matrix).
#
# To add new data structures, each of the below functions must be updated.

must.be.supported.data.type <- function(X) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(X)), ".")
  if (!(is.matrix(X)
        || inherits(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || (inherits(X, "lowrank") || is.udv(X))
        || (inherits(X, "lrps") || is.UVS(X))))
    stop(error.msg)
}

get.data.dims <- function(X) {
  if (inherits(X, "lowrank")) {
    return(sapply(X, nrow))
  } else if (inherits(X, "lrps")) {
    return(dim(X[["sparse"]]))
  } else {
    return(dim(X))
  }
}

get.data.dimnames <- function(X) {
  if (inherits(X, "lowrank")) {
    return(lapply(X, rownames))
  } else if (inherits(X, "lrps")) {
    return(lapply(X[["lowrank"]], rownames))
  } else {
    return(dimnames(X))
  }
}

as.fullrank <- function(X) {
  if (inherits(X, "lowrank")) {
    return(lowrank.expand(X))
  } else if (inherits(X, "lrps")) {
    return(lowrank.expand(X[["lowrank"]]) + X[["sparse"]])
  } else {
    return(X)
  }
}

nmode.ops.error <- paste("N-mode products are not yet implemented for",
                         "objects of that class and/or dimension.")

# X is a m_1 x m_2 (x m_3) matrix (array) and v is an m_n-vector. The n-mode
#   product multiplies the jth n-slice of X by v_j and then sums the slices
#   together. For example, the 3-mode product of a 3-dimensional tensor is
#   \sum_k v_k X_{ijk}. v can be set to 1 to indicate an m_n-vector of all
#   ones.
#
#' @importFrom Matrix colSums rowSums t crossprod
#'
nmode.prod.vec <- function(X, v, n) {
  # Matrices and sparse matrices:
  if (is.matrix(X) || inherits(X, "Matrix")) {
    if (n == 1 && identical(v, 1))
      return(colSums(X))
    if (n == 1)
      return(as.vector(v %*% X))
    if (n == 2 && identical(v, 1))
      return(rowSums(X))
    if (n == 2)
      return(as.vector(X %*% v))
  }

  # Low-rank matrices:
  if (inherits(X, "lowrank") && length(X) == 2) {
    if (identical(v, 1))
      u <- colSums(X[[n]])
    else
      u <- t(v %*% X[[n]])

    return(as.vector(X[[-n]] %*% u))
  }

  # Low-rank plus sparse:
  if (inherits(X, "lrps")) {
    return(
      nmode.prod.vec(X[["lowrank"]], v, n) +
        nmode.prod.vec(X[["sparse"]], v, n)
    )
  }

  # 3-dimensional tensors:
  if (is.array(X) && length(dim(X) == 3)) {
    if (n == 1 && identical(v, 1))
      return(apply(X, 3, colSums))
    if (n == 1)
      return(apply(X, 3, FUN = function(M) {t(v) %*% M}))
    if (n == 2 && identical(v, 1))
      return(apply(X, 3, rowSums))
    if (n == 2)
      return(apply(X, 3, FUN = function(M) {M %*% v}))
    if (n == 3 && identical(v, 1))
      return(apply(X, 2, rowSums))
    if (n == 3)
      return(apply(X, 2, FUN = function(M) {M %*% v}))
  }

  # Low-rank tensors (returns a low-rank matrix):
  if (inherits(X, "lowrank") && length(X) == 3) {
    if (identical(v, 1))
      u <- colSums(X[[n]])
    else
      u <- as.vector(v %*% X[[n]])

    retval <- X[-n]
    smaller.dim <- which.min(c(nrow(retval[[1]]), nrow(retval[[2]])))
    retval[[smaller.dim]] <- t(t(retval[[smaller.dim]]) * u)
    class(retval) <- "lowrank"
    return(retval)
  }

  stop(nmode.ops.error)
}

# X is a m_1 x m_2 (x m_3) matrix (array), r1 is a r1 object whose dimension
#   is one less than X, and n indicates the dimension that is omitted. For
#   matrices, this is simply the usual matrix product of X and the vector
#   in r1. For tensors, it is the n-mode product of X and the rank-one matrix
#   formed by the outer product of the vectors in r1. For example, if
#   r1 = list(u, v) and n = 1, then \sum_{j, k} u_j v_k X_{ijk} is returned.
#   Any of the list elements in r1 can be set to 1 to indicate a vector of all
#   ones.
#
nmode.prod.r1 <- function(X, r1, n) {
  if (is.null(X))
    return(0)

  if (identical(X, 1)) {
    return(r1.sum(r1))
  }

  if (is.matrix(X)
      || inherits(X, "Matrix")
      || (inherits(X, "lowrank") && length(X) == 2)
      || inherits(X, "lrps")) {
    return(nmode.prod.vec(X, unlist(r1), (1:2)[-n]))
  }

  if ((is.array(X) && (length(dim(X)) == 3))
      || (inherits(X, "lowrank") && length(X) == 3)) {
    ns <- (1:3)[-n]
    # Go backwards (otherwise indices get messed up).
    return(nmode.prod.vec(nmode.prod.vec(X, r1[[2]], ns[2]), r1[[1]], ns[1]))
  }

  stop(nmode.ops.error)
}

nmode.prod.lowrank <- function(X, Y, n) {
  # Matrices and sparse matrices:
  if (is.matrix(X) || inherits(X, "Matrix")) {
    if (n == 1)
      return(as.matrix(X %*% Y[[1]]))
    if (n == 2)
      return(t(as.matrix(t(Y[[1]]) %*% X)))
  }

  # Low-rank matrix representations:
  if (inherits(X, "lowrank") && (length(X) == 2)) {
    if (n == 1)
      return(X[[1]] %*% (t(X[[2]]) %*% Y[[1]]))
    if (n == 2)
      return(X[[2]] %*% (t(X[[1]]) %*% Y[[1]]))
  }

  # Low-rank plus sparse:
  if (inherits(X, "lrps")) {
    return(nmode.prod.lowrank(X[["lowrank"]], Y, n) +
             nmode.prod.lowrank(X[["sparse"]], Y, n))
  }

  if ((is.array(X) && (length(dim(X)) == 3))
      || (inherits(X, "lowrank") && (length(X) == 3))) {
    stop("N-mode products of tensors with matrices have not yet been implemented.")
  }

  stop(nmode.ops.error)
}

# The following calculates the n-mode product between the matrix or tensor
#   X * LF (where LF is low-rank and multiplication is elementwise) and an r1
#   object. The idea is both to avoid forming a large matrix (or tensor) and
#   to be able to take advantage of any sparsity in X (since LF will not in
#   general be sparse, simply performing elementwise multiplication will
#   destroy sparsity in X).
#
premult.nmode.prod.r1 <- function(Z, X, r1, n) {
  if (is.null(X))
    return(0)

  if (inherits(X, "lowrank"))
    return(premult.lowrank.nmode.prod.r1(Z, X, r1, n))

  if (inherits(Z, "lowrank"))
    return(premult.lowrank.nmode.prod.r1(X, Z, r1, n))

  if (inherits(X, "lrps")) {
    return(premult.nmode.prod.r1(Z, X[["lowrank"]], r1, n) +
             premult.nmode.prod.r1(Z, X[["sparse"]], r1, n))
  }

  if (inherits(Z, "lrps")) {
    return(premult.nmode.prod.r1(Z[["lowrank"]], X, r1, n) +
             premult.nmode.prod.r1(Z[["sparse"]], X, r1, n))
  }

  return(nmode.prod.r1(Z * X, r1, n))
}

premult.lowrank.nmode.prod.r1 <- function(Z, lowrank, r1, n) {
  if (identical(Z, 1)) {
    return(nmode.prod.r1(lowrank, r1, n))
  }

  # Matrices:
  if (length(lowrank) == 2) {
    u <- unlist(r1) * lowrank[[-n]]

    if (inherits(Z, "lowrank") && n == 1)
      return(rowSums(lowrank[[1]] * (Z[[1]] %*% crossprod(Z[[2]], u))))
    if (inherits(Z, "lowrank") && n == 2)
      return(rowSums(lowrank[[2]] * (Z[[2]] %*% t(crossprod(u, Z[[1]])))))
    if (n == 1)
      return(rowSums(lowrank[[1]] * (Z %*% u)))
    if (n == 2)
      return(colSums(t(lowrank[[2]]) * (t(u) %*% Z)))
  }

  # 3-dimensional tensors:
  if (length(lowrank) == 3) {
    ns <- (1:3)[-n]
    # Sum over third (or second) index. Get k m_1 x m_2 (or m_1 x m_3) matrices
    #   stacked on top of one another.
    tmp <- nmode.prod.vec(Z, r1[[2]] * lowrank[[ns[2]]], ns[[2]])
    if (inherits(tmp, "lowrank"))
      tmp <- lowrank.expand(tmp)
    # Elementwise multiply by an m_1 x k or m_2 x k matrix, then sum over the
    #   m_ index.
    if (ns[[1]] == 1) {
      tmp <- tmp * as.vector(lowrank[[1]] * r1[[1]])
      tmp <- colSums(matrix(tmp, nrow = dim(Z)[1]))
      # Have m_n k vectors concatenated one after another.
      return(rowSums(lowrank[[n]]
                             * matrix(tmp, nrow = dim(Z)[n], byrow = TRUE)))
    } else { # ns[[1]] == 2
      tmp <- tmp * rep(t(lowrank[[2]] * r1[[1]]), each = dim(Z)[1])
      tmp <- rowSums(tmp)
      # Have k m_n vectors concatenated one after another.
      return(rowSums(lowrank[[n]] * matrix(tmp, nrow = dim(Z)[n])))
    }
  }

  stop(nmode.ops.error)
}

is.udv <- function(X) {
  if (!is.list(X))
    return(FALSE)

  # Must have fields d, u, and v.
  if (is.null(X$d) || is.null(X$u) || is.null(X$v))
    return(FALSE)

  # Check u and v.
  if (!(is.matrix(X$u) || inherits(X$u, "Matrix")) ||
      !(is.matrix(X$v) || inherits(X$v, "Matrix")))
    return(FALSE)
  if (!identical(ncol(X$u), ncol(X$v)))
    return(FALSE)
  if (anyNA(X$u) || anyNA(X$v))
    return(FALSE)

  # Check d.
  if (!is.numeric(X$d) || (length(X$d) < ncol(X$u)))
    return(FALSE)

  return(TRUE)
}

is.UVS <- function(X) {
  if (!is.list(X))
    return(FALSE)

  # Must have fields U, V, and S.
  if (is.null(X$U) || is.null(X$V) || is.null(X$S))
    return(FALSE)

  # Check U and V.
  if (!(is.matrix(X$U) || inherits(X$U, "Matrix")) ||
      !(is.matrix(X$V) || inherits(X$V, "Matrix")))
    return(FALSE)
  if (!identical(ncol(X$U), ncol(X$V)))
    return(FALSE)
  if (anyNA(X$U) || anyNA(X$v))
    return(FALSE)

  # Check S.
  if (!(is.matrix(X$S) || inherits(X$S, "Matrix")))
    return(FALSE)
  if (!identical(nrow(X$U), nrow(X$S)))
    return(FALSE)
  if (!identical(nrow(X$V), ncol(X$S)))
    return(FALSE)

  return(TRUE)
}
