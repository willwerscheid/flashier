must.be.supported.data.type <- function(X,
                                        allow.null = TRUE,
                                        allow.scalar = FALSE,
                                        allow.vector = FALSE) {
  if (!(is(X, "flash.data")
        || is.matrix(X)
        || is(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || (allow.null && is.null(X))
        || (allow.scalar && is.vector(X) && length(X) == 1)
        || (allow.vector && is.vector(X))))
    stop()
}

must.be.integer <- function(x, lower = NULL, upper = NULL, allow.null = TRUE) {
  if (is.null(x)) {
    if (!allow.null)
      stop()
  } else if (!(is.numeric(x)
               && is.finite(x)
               && as.integer(x) == x
               && (is.null(lower) || x >= lower)
               && (is.null(upper) || x <= upper)))
    stop()
}

must.be.named.list <- function(x) {
  if (!is.null(x) && (!is.list(x) || is.null(names(x))))
    stop()
}

dims.must.match <- function(X, Y, n = NULL) {
  if (is.null(n)) {
    if (!is.null(X) && !is.null(Y) && !identical(dim(X), dim(Y)))
      stop()
  } else {
    if (n == 0 && length(Y) != 1)
      stop()
    if (n > 0 && (!is.vector(Y) || !identical(length(Y), dim(X)[n])))
      stop()
  }
}
