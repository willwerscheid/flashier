#   Similar to r1 objects, lowrank objects are lists of matrices whose columns
# each correspond to a single rank-one matrix or tensor. For example,
# lowrank = list(U, V, W) describes a tensor whose ijk-entry is equal to
# \sum_{\ell} U_{i \ell} V_{j \ell} W_{k \ell}. A lowrank class is defined but
# these operations do not usually check the class of their arguments.

as.lowrank <- function(r1) {
  if (is.null(r1))
    return(NULL)
  if (inherits(r1, "r1")) {
    lowrank <- lapply(r1, matrix, ncol = 1)
    class(lowrank) <- "lowrank"
    return(lowrank)
  }
  stop("as.lowrank not defined for that class.")
}

# This function is slow and should only be used when calculating residuals.
lowrank.expand <- function(lowrank) {
  if (is.null(lowrank))
    return(0)

  if (length(lowrank) == 2)
    return(tcrossprod(lowrank[[1]], lowrank[[2]]))

  # TODO: This is probably not the best way to do this...
  K   <- ncol(lowrank[[1]])
  r1s <- lapply(lowrank, function(M) {split(M, col(M))})
  r1s <- lapply(1:K, function(k) {lapply(r1s, function(r1) {r1[[k]]})})

  res <- r1.expand(r1s[[1]])
  if (K > 1)
    for (k in 2:K)
      res <- res + r1.expand(r1s[[k]])
  return(res)
}

lowrank.subset <- function(lowrank, n, subset) {
  if (is.null(lowrank))
    return(NULL)
  lowrank[[n]] <- lowrank[[n]][subset, , drop = FALSE]
  return(lowrank)
}

lowrank.drop.k <- function(lowrank, k) {
  if (is.null(lowrank) || all(1:ncol(lowrank[[1]]) %in% k))
    return(NULL)
  lowrank <- lapply(lowrank, function(X) X[, -k, drop = FALSE])
  class(lowrank) <- "lowrank"
  return(lowrank)
}

lowrank.sc.mult <- function(lowrank, x) {
  if (is.null(lowrank))
    return(NULL)
  lowrank[[1]] <- x * lowrank[[1]]
  return(lowrank)
}

# This squares each entry of L and each entry of F. It does not square the
#   entries of LF'.
lowrank.square <- function(lowrank) {
  if (is.null(lowrank))
    return(NULL)
  lowrank2 <- lapply(lowrank, function(M) {M^2})
  class(lowrank2) <- "lowrank"
  return(lowrank2)
}

lowranks.prod <- function(lr1, lr2, broadcast = FALSE) {
  if (is.null(lr1) || is.null(lr2))
    return(NULL)
  if (broadcast)
    lr1 <- lapply(lr1, as.vector)
  product <- mapply(`*`, lr1, lr2, SIMPLIFY = FALSE)
  class(product) <- "lowrank"
  return(product)
}

lowranks.combine <- function(lr1, lr2) {
  if (is.null(lr1))
    return(lr2)
  if (is.null(lr2))
    return(lr1)
  lowrank <- mapply(cbind, lr1, lr2, SIMPLIFY = FALSE)
  class(lowrank) <- "lowrank"
  return(lowrank)
}

lowrank.delta.mat <- function(new.lr, old.lr) {
  k <- ncol(new.lr[[1]])
  lowrank <- mapply(cbind, new.lr, old.lr, SIMPLIFY = FALSE)
  if (!is.null(k) && k > 0)
    lowrank[[1]][, (k + 1):(2 * k)] <- -lowrank[[1]][, (k + 1):(2 * k)]
  class(lowrank) <- "lowrank"
  return(lowrank)
}

elemwise.prod.lowrank.r1 <- function(lowrank, r1) {
  lr.prod <- mapply(`*`, lowrank, r1, SIMPLIFY = FALSE)
  class(lr.prod) <- "lowrank"
  return(lr.prod)
}
