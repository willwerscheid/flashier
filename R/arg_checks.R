# TODO: error messages

must.be.integer <- function(x, lower = NULL, upper = NULL, allow.null = TRUE) {
  if (is.null(x) && !allow.null)
    stop()
  if (!is.null(x) && (!(is.numeric(x)
                        && is.finite(x)
                        && as.integer(x) == x
                        && (is.null(lower) || x >= lower)
                        && (is.null(upper) || x <= upper))))
    stop()
}

must.be.named.list <- function(x) {
  if (!is.null(x) && (!is.list(x) || is.null(names(x))))
    stop()
}

must.be.supported.data.type <- function(X,
                                        allow.null = TRUE,
                                        allow.vector = FALSE) {
  if (!(is(X, "flash.data")
        || is.matrix(X)
        || is(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || (allow.null && is.null(X))
        || (allow.vector && is.vector(X))))
    stop()
}

must.be.valid.var.type <- function(x, data.dim, allow.null = TRUE) {
  if (is.null(x) && !allow.null)
    stop()
  for (val in x) {must.be.integer(val, lower = 0, upper = data.dim)}
  # No repeats.
  if (!identical(x, unique(x)))
    stop()
  # If zero appears (meaning "constant" var.type), it must be alone.
  if ((0 %in% x) && (length(x) > 1))
    stop()
}

dims.must.match <- function(X, Y, Y.dim = NULL) {
  # If Y.dim is NULL, then Y must be a matrix or array.
  if (is.null(Y.dim)) {
    if (!is.null(X) && !is.null(Y) && !identical(dim(X), dim(Y)))
      stop()
  } else {
    # If Y.dim is zero, then Y must be a scalar.
    if (Y.dim == 0 && (!is.vector(Y) || (length(Y) != 1)))
      stop()
    # Otherwise, Y must be a vector that can be aligned with X.
    if (Y.dim > 0 && (!is.vector(Y) || (length(Y) != dim(X)[Y.dim])))
      stop()
  }
}
