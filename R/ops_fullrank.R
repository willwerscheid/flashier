# Fullrank operations are implemented for the same types of objects as
#   n-mode products.

fullrank.ops.error <- function(fn.name, object) {
  paste(fn.name, "is not yet implemented for objects of the same class and/or",
        "dimension as object", deparse(substitute(object)))
}

fullrank.subset <- function(X, n, subset) {
  if(is.null(X))
    return(NULL)

  if(identical(X, 1))
    return(1)

  if (is.matrix(X) || is(X, "Matrix")) {
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

  stop(fullrank.ops.error("fullrank.subset", X))
}
