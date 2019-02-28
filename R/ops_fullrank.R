#   Fullrank operations are implemented for the same types of objects as
# n-mode products.

fullrank.ops.error <- function(fn.name, object) {
  paste(fn.name, "is not yet implemented for objects of that class and/or",
        "dimension.")
}

full.or.lowrank.subset <- function(X, n, subset) {
  if (inherits(X, "lowrank"))
    return(lowrank.subset(X, n, subset))
  return(fullrank.subset(X, n, subset))
}

fullrank.subset <- function(X, n, subset) {
  if(is.null(X))
    return(NULL)

  if(identical(X, 1))
    return(1)

  if (is.matrix(X) || inherits(X, "Matrix")) {
    if (n == 1)
      return(X[subset, , drop = FALSE])
    if (n == 2)
      return(X[, subset, drop = FALSE])
  }

  if (is.array(X) && length(dim(X) == 3)) {
    if (n == 1)
      return(X[subset, , , drop = FALSE])
    if (n == 2)
      return(X[, subset, , drop = FALSE])
    if (n == 3)
      return(X[, , subset, drop = FALSE])
  }

  stop(fullrank.ops.error("fullrank.subset"))
}

elemwise.prod.fullrank.r1 <- function(X, r1) {
  return(X * r1.expand(r1))
}

elemwise.prod.fullrank.lowrank <- function(X, lr) {
  return(X * lowrank.expand(lr))
}
