#   r1 objects are simply lists of vectors. For example, r1 = list(u, v, w)
# describes a tensor with ijk-entry equal to u_i v_j w_k. An r1 class is
# defined but these operations do not check the class of their arguments.

r1.ops.error <- function(fn.name, dim) {
  paste(fn.name, "is not yet implemented for r1 objects of dimension", dim)
}

r1.expand <- function(r1) {
  if (length(r1) == 1)
    return(r1[[1]])
  if (length(r1) == 2)
    return(outer(r1[[1]], r1[[2]]))
  if (length(r1) == 3)
    return(r1[[1]] %o% r1[[2]] %o% r1[[3]])
  stop(r1.ops.error("r1.expand", length(r1)))
}

r1.subset <- function(r1, n, subset) {
  r1[[n]] <- r1[[n]][subset]
  return(r1)
}

r1.drop.dim <- function(r1, n) {
  r1 <- r1[-n]
  class(r1) <- "r1"
  return(r1)
}

r1.square <- function(r1) {
  r12 <- lapply(r1, function(v) {v^2})
  class(r12) <- "r1"
  return(r12)
}

r1.sum <- function(r1) {
  return(prod(sapply(r1, sum)))
}

r1.zeros <- function(flash) {
  r1 <- lapply(get.dims(flash), function(dim) rep(0, dim))
  class(r1) <- "r1"
  return(r1)
}

r1.ones <- function(flash) {
  r1 <- as.list(rep(1, get.dim(flash) - 1))
  class(r1) <- "r1"
  return(r1)
}

#' @importFrom stats rnorm
r1.random <- function(dims, dim.signs = NULL) {
  r1 <- lapply(dims, function(dim) {rnorm(dim) / sqrt(dim)})
  for (n in which(dim.signs == 1)) {
    r1[[n]] <- abs(r1[[n]])
  }
  for (n in which(dim.signs == -1)) {
    r1[[n]] <- -abs(r1[[n]])
  }
  class(r1) <- "r1"
  return(r1)
}

r1.to.lowrank <- function(r1, flash) {
  lowrank <- mapply(r1, get.dims(flash),
                    FUN = function(vals, dim) {
                      matrix(vals, ncol = 1, nrow = dim)
                    }, SIMPLIFY = FALSE)
  class(lowrank) <- "lowrank"
  return(lowrank)
}
