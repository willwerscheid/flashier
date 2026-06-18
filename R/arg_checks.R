must.be.numeric <- function(x, allow.infinite = TRUE, allow.null = TRUE) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(x)), ".")
  if (is.null(x) && !allow.null)
    stop(error.msg)
  if (!is.null(x) && is.infinite(x) && !allow.infinite)
    stop(error.msg)
  if (!is.null(x) && (!(is.numeric(x) && length(x) == 1)))
    stop(error.msg)
}

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

must.be.valid.kset <- function(flash, kset) {
  for (k in kset) {
    must.be.integer(k, lower = 1, upper = get.n.factors(flash),
                    allow.null = FALSE)
  }
  if (!identical(kset, unique(kset))) {
    stop("kset must not contain repeated elements.")
  }
}

must.be.list.of.named.lists <- function(x, valid.fields) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(x)), ".")
  if (!is.list(x) || !all(sapply(x, is.list)))
    stop(error.msg)
  if (!all(unlist(lapply(x, names)) %in% valid.fields))
    stop(error.msg)
}

must.be.supported.S.type <- function(X) {
  error.msg <- paste0("Invalid argument to ", deparse(substitute(X)), ".")
  if (!(is.matrix(X)
        || inherits(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || is.null(X)
        || is.vector(X))) {
    stop(error.msg)
  }
}

must.be.compatible.data.types <- function(X, Y) {
  error.msg <- paste0("If either ", deparse(substitute(X)), " or ",
                      deparse(substitute(Y)), " is of class Matrix, then both",
                      " must be.")
  if ((inherits(X, "Matrix") && is.matrix(Y))
      || (is.matrix(X) && inherits(Y, "Matrix")))
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

must.not.have.zero.slices <- function(Y, S, var.type) {
  # Skip this test for non-matrix data structures.
  if (!(is.matrix(Y) || inherits(Y, "Matrix")))
    return()

  # The purpose of the check is to prevent precision estimates from going to
  #   zero, so we can skip it if S > 0 or if var.type = 0.
  if (!is.null(S) && all(S > 0))
    return()
  if (!is.null(var.type) && length(var.type) == 1 && var.type == 0)
    return()

  error.msg <- paste("The data matrix must not have any rows or",
                     "columns whose entries are either identically zero",
                     "or all missing.")

  nz <- (Y != 0)
  for (n in 1:length(dim(Y))) {
    n.nonzero <- nmode.prod.vec(nz, 1, n)
    if (any(n.nonzero == 0) && (is.null(var.type) || n %in% var.type))
      stop(error.msg)
  }
}

dims.must.match <- function(X, Y, Y.dim = NULL) {
  error.msg <- paste("Dimensions of", deparse(substitute(X)), "and",
                     deparse(substitute(Y)), "do not match.")
  if (inherits(X, "lowrank")) {
    dim.X <- sapply(X, nrow)
  } else {
    dim.X <- dim(X)
  }
  # If Y.dim is NULL, then Y must be a matrix or array.
  if (is.null(Y.dim)) {
    if (!is.null(X) && !is.null(Y) && !identical(dim.X, get.data.dims(Y)))
      stop(error.msg)
  } else {
    # If Y.dim is zero, then Y must be a scalar.
    if (Y.dim == 0 && (!is.vector(Y) || (length(Y) != 1)))
      stop(error.msg)
    # Otherwise, Y must be a vector that can be aligned with X.
    if (Y.dim > 0 && (!is.vector(Y) || (length(Y) != dim.X[Y.dim])))
      stop(error.msg)
  }
}
