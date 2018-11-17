# Basic operations. These work for matrices, sparse matrices, and
#   3-dimensional tensors (stored as arrays). I would also like to implement
#   them for sparse tensors and tensors of arbitrary dimension.

basic.ops.error <- paste("Basic operations not yet implemented for objects",
                         "of this class and dimension.")

# N-mode products -------------------------------------------------------------

# X is a m_1 x m_2 (x m_3) matrix (array) and v is an m_n-vector. The n-mode
#   product multiplies the jth n-slice of X by v_j and then sums the slices
#   together. For example, the 3-mode product of a 3-dimensional tensor is
#   \sum_k v_k X_{ijk}.
# v can be set to 1 to indicate an m_n-vector of all ones.
#
nmode.prod.vec <- function(X, v, n) {
  # Matrices and sparse matrices:
  if (is.matrix(X) || is(X, "Matrix")) {
    if (n == 1 && identical(v, 1))
      return(Matrix::colSums(X))
    if (n == 1)
      return(as.vector(v %*% X))
    if (n == 2 && identical(v, 1))
      return(Matrix::rowSums(X))
    if (n == 2)
      return(as.vector(X %*% v))
  }

  # 3-dimensional tensors:
  if (is.array(X) && length(dim(X) == 3)) {
    if (n == 1 && identical(v, 1))
      return(apply(X, 3, Matrix::colSums))
    if (n == 1)
      return(apply(X, 3, FUN = function(M) {t(v) %*% M}))
    if (n == 2 && identical(v, 1))
      return(apply(X, 3, Matrix::rowSums))
    if (n == 2)
      return(apply(X, 3, FUN = function(M) {M %*% v}))
    if (n == 3 && identical(v, 1))
      return(apply(X, 2, Matrix::rowSums))
    if (n == 3)
      return(apply(X, 2, FUN = function(M) {M %*% v}))
  }

  stop(basic.ops.error)
}

# X is a m_1 x m_2 (x m_3) matrix (array), v is a r1 object whose dimension
#   is one less than X, and n indicates the dimension that is omitted. For
#   matrices, this is simply the usual matrix product of X and the vector
#   in r1. For tensors, it is the n-mode product of X and the rank-one matrix
#   formed by the outer product of the vectors in r1. For example, if
#   r1 = list(u, v) and n = 1, then \sum_{j, k} u_j v_k X_{ijk} is returned.
# Any of the list elements in r1 can be set to 1 to indicate a vector of all
#   ones.
#
fullrank.nmode.prod.r1 <- function(X, r1, n) {
  if (identical(X, 1)) {
    return(r1.sum(r1))
  }

  if (is.matrix(X) || is(X, "Matrix")) {
    return(nmode.prod.vec(X, unlist(r1), (1:2)[-n]))
  }

  if (is.array(X) && length(dim(X) == 3)) {
    ns <- (1:3)[-n]
    # Go backwards (otherwise indices get messed up):
    return(nmode.prod.vec(nmode.prod.vec(X, r1[[2]], ns[2]), r1[[1]], ns[1]))
  }

  stop(basic.ops.error)
}

# The following function has the same purpose as nmode.prod.r1, except that
#   it takes a lowrank object as an argument rather than a full matrix (or
#   tensor). The idea is to avoid computing with a large m_1 x m_2 matrix (or
#   m_1 x m_2 x m_3 tensor).
#
lowrank.nmode.prod.r1 <- function(lowrank, r1, n) {
  # Matrices:
  if (length(lowrank) == 2)
    return(as.vector(lowrank[[n]]
                     %*% nmode.prod.vec(lowrank[[(1:2)[-n]]], r1[[1]], 1)))

  # 3-dimensional tensors:
  if (length(lowrank) == 3) {
    ns <- (1:3)[-n]
    return(as.vector(lowrank[[n]]
                     %*% (nmode.prod.vec(lowrank[[ns[1]]], r1[[1]], 1)
                          * nmode.prod.vec(lowrank[[ns[2]]], r1[[2]], 1))))
  }

  stop(basic.ops.error)
}

nmode.prod.r1 <- function(X, r1, n) {
  if (is(X, "lowrank"))
    return(lowrank.nmode.prod.r1(X, r1, n))
  return(fullrank.nmode.prod.r1(X, r1, n))
}

# This calculates the n-mode product between the matrix or tensor X * LF
#   (where LF is low-rank and multiplication is elementwise) and an r1
#   object. The idea is both to avoid forming a large matrix (or tensor) and
#   to be able to take advantage of any sparsity in X (since LF will not in
#   general be sparse, simply performing elementwise multiplication will
#   destroy sparsity in X).
#
premult.lowrank.nmode.prod.r1 <- function(Z, lowrank, r1, n) {
  if (identical(Z, 1)) {
    return(lowrank.nmode.prod.r1(lowrank, r1, n))
  }

  # Matrices and sparse matrices:
  if (length(lowrank) == 2) {
    if (n == 1)
      return(Matrix::rowSums(lowrank[[1]]
                             * (Z %*% (unlist(r1) * lowrank[[2]]))))
    if (n == 2)
      return(Matrix::colSums(t(lowrank[[2]])
                             * (t(unlist(r1) * lowrank[[1]]) %*% Z)))
  }

  if (length(lowrank) == 3) {
    ns <- (1:3)[-n]
    # Sum over third (or second) index. Get k m_1 x m_2 (or m_1 x m_3) matrices
    #   stacked on top of one another:
    tmp <- nmode.prod.vec(Z, r1[[2]] * lowrank[[ns[2]]], ns[[2]])
    # Elementwise multiply by an m_1 x k or m_2 x k matrix, then sum over the
    #   m_ index.
    if (ns[[1]] == 1) {
      tmp <- tmp * as.vector(lowrank[[1]] * r1[[1]])
      tmp <- Matrix::colSums(matrix(tmp, nrow = dim(Z)[1]))
      # Have m_n k vectors concatenated one after another:
      return(Matrix::rowSums(lowrank[[n]]
                             * matrix(tmp, nrow = dim(Z)[n], byrow = TRUE)))
    } else { # ns[[1]] == 2
      tmp <- tmp * rep(t(lowrank[[2]] * r1[[1]]), each = dim(Z)[1])
      tmp <- Matrix::rowSums(tmp)
      # Have k m_n vectors concatenated one after another:
      return(Matrix::rowSums(lowrank[[n]] * matrix(tmp, nrow = dim(Z)[n])))
    }
  }

  stop(basic.ops.error)
}

premult.nmode.prod.r1 <- function(Z, X, r1, n) {
  if (is.null(X))
    return(0)

  if (is(X, "lowrank"))
    return(premult.lowrank.nmode.prod.r1(Z, X, r1, n))

  return(fullrank.nmode.prod.r1(Z * X, r1, n))
}


# Functions on r1 objects -----------------------------------------------------

# r1 objects are simply lists of vectors. For example, r1 = list(u, v, w)
#   describes a tensor with ijk-entry equal to u_i v_j w_k.
#
# TODO rewrite this one
r1.expand <- function(r1) {
  tmp <- r1[[1]]
  for (i in 2:length(r1)) {
    tmp <- tmp %o% r1[[i]]
  }
  return(tmp)
}

r1.sum <- function(r1) {
  return(prod(sapply(r1, sum)))
}

r1.square <- function(r1) {
  r12 <- lapply(r1, function(v) {v^2})
  class(r12) <- "r1"
  return(r12)
}

r1.ones <- function(flash) {
  r1 <- as.list(rep(1, get.dim(flash) - 1))
  class(r1) <- "r1"
  return(r1)
}

r1.random <- function(dim) {
  r1 <- lapply(dim, function(d) {rnorm(d) / sqrt(d)})
  class(r1) <- "r1"
  return(r1)
}

# Functions on lowrank objects ------------------------------------------------

# Similar to r1 objects, lowrank objects are lists of matrices whose columns
#   each correspond to a single rank-one matrix or tensor. For example,
#   lowrank = list(U, V, W) describes a tensor whose ijk-entry is equal to
#   \sum_{\ell} U_{i \ell} V_{j \ell} W_{k \ell}.

# TODO: maybe rewrite, but this is only used rarely.
lowrank.expand <- function(lowrank) {
  if (is.null(lowrank))
    return(0)
  if (is.matrix(lowrank) || is(lowrank, "Matrix"))
    return(tcrossprod(lowrank[[1]], lowrank[[2]]))

  r1s <- list()
  for (k in 1:ncol(lowrank[[1]])) {
    r1s[[k]] <- lapply(lowrank, function(M) {M[, k]})
  }
  res <- r1.expand(r1s[[1]])
  if (length(r1s) > 1)
    for (k in 2:length(r1s))
      res <- res + r1.expand(r1s[[k]])
  return(res)
}

lowrank.square <- function(lowrank) {
  if (is.null(lowrank))
    return(NULL)
  lowrank2 <- lapply(lowrank, function(M) {M^2})
  class(lowrank2) <- "lowrank"
  return(lowrank2)
}

lowrank.sc.mult <- function(lowrank, x) {
  if (is.null(lowrank))
    return(0)
  lowrank[[1]] <- x * lowrank[[1]]
  return(lowrank)
}

lowrank.drop.k <- function(lowrank, k) {
  if (is.null(lowrank))
    return(NULL)
  lowrank <- lapply(lowrank, function(X) X[, -k, drop = FALSE])
  class(lowrank) <- "lowrank"
  return(lowrank)
}

as.lowrank <- function(r1) {
  if (is.null(r1))
    return(NULL)
  if (is(r1, "r1")) {
    lowrank <- lapply(r1, function(n) {matrix(n, ncol = 1)})
    class(lowrank) <- "lowrank"
    return(lowrank)
  } else {
    stop("as.lowrank not defined for that class.")
  }
}

lowranks.prod <- function(lr1, lr2, broadcast = FALSE) {
  if (is.null(lr1) || is.null(lr2))
    return(NULL)
  if (broadcast)
    lr1 <- lapply(lr1, as.vector)
  prod <- mapply(`*`, lr1, lr2)
  class(prod) <- "lowrank"
  return(prod)
}

lowranks.combine <- function(lr1, lr2) {
  if (is.null(lr1))
    return(lr2)
  if (is.null(lr2))
    return(lr1)
  lowrank <- mapply(cbind, lr1, lr2)
  class(lowrank) <- "lowrank"
  return(lowrank)
}

elemwise.prod.lowrank.r1 <- function(lowrank, r1) {
  lr.prod <- mapply(lowrank, r1, FUN = `*`)
  class(lr.prod) <- "lowrank"
  return(lr.prod)
}

lowrank.delta.mat <- function(new.lr, old.lr) {
  n <- ncol(new.lr[[1]])
  lowrank <- mapply(cbind, new.lr, old.lr)
  lowrank[[1]][, (n + 1):(2 * n)] <- -lowrank[[1]][, (n + 1):(2 * n)]
  class(lowrank) <- "lowrank"
  return(lowrank)
}
