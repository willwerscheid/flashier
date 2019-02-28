must.be.integer <- function(x, lower = NULL, upper = NULL, allow.null = TRUE) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(x)), ".")
  if (is.null(x) && !allow.null)
    stop(error.msg)
  if (!is.null(x) && (!(is.numeric(x)
                        && is.finite(x)
                        && as.integer(x) == x
                        && (is.null(lower) || x >= lower)
                        && (is.null(upper) || x <= upper))))
    stop(error.msg)
}

must.be.named.list <- function(x) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(x)), ".")
  if (!is.list(x) || (length(x) > 0 && is.null(names(x))))
    stop(error.msg)
}

must.be.supported.data.type <- function(X,
                                        allow.null = TRUE,
                                        allow.vector = FALSE) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(X)), ".")
  if (!(inherits(X, "flash.data")
        || is.matrix(X)
        || inherits(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || (allow.null && is.null(X))
        || (allow.vector && is.vector(X))))
    stop(error.msg)
}

must.be.valid.var.type <- function(x, data.dim, allow.null = TRUE) {
  error.msg <- "Invalid var.type."
  if (is.null(x) && !allow.null)
    stop(error.msg)
  for (val in x) {must.be.integer(val, lower = 0, upper = data.dim)}
  # No repeats.
  if (!identical(x, unique(x)))
    stop(error.msg)
  # If zero appears (meaning "constant" var.type), it must be alone.
  if ((0 %in% x) && (length(x) > 1))
    stop(error.msg)
}

must.not.have.zero.slices <- function(Y) {
  nz <- (Y != 0)
  for (n in 1:length(dim(Y))) {
    n.nonzero <- nmode.prod.vec(nz, 1, n)
    if (any(n.nonzero == 0))
      stop("The data matrix must not have any slices (rows, columns) whose",
           " entries are identically zero.")
  }
}

dims.must.match <- function(X, Y, Y.dim = NULL) {
  error.msg <- paste("Dimensions of", deparse(substitute(X)), "and",
                     deparse(substitute(Y)), "do not match.")
  # If Y.dim is NULL, then Y must be a matrix or array.
  if (is.null(Y.dim)) {
    if (!is.null(X) && !is.null(Y) && !identical(dim(X), dim(Y)))
      stop(error.msg)
  } else {
    # If Y.dim is zero, then Y must be a scalar.
    if (Y.dim == 0 && (!is.vector(Y) || (length(Y) != 1)))
      stop(error.msg)
    # Otherwise, Y must be a vector that can be aligned with X.
    if (Y.dim > 0 && (!is.vector(Y) || (length(Y) != dim(X)[Y.dim])))
      stop(error.msg)
  }
}
