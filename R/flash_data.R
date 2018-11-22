set.flash.data <- function(data, S = NULL, tau = NULL, tau.dim = NULL) {
  if (is(data, "flash.data")) {
    for (arg in list(S, tau, tau.dim)) {
      warning("Ignoring ", deparse(substitute(arg)))
    }
  } else {
    must.not.supply.both(S, tau)
    must.not.supply.both(S, tau.dim)
    must.supply.neither.or.both(tau, tau.dim)
    must.be.supported.data.type(data, allow.null = FALSE)
    must.be.supported.data.type(S, allow.scalar = TRUE)
    must.be.supported.data.type(tau, allow.vector = TRUE)
    must.be.valid.integer(tau.dim, lower = 0, upper = length(dim(data)))

    flash.data        <- list()
    class(flash.data) <- "flash.data"
    flash.data$Y      <- data
    if (anyNA(data)) {
      flash.data$Y[is.na(data)] <- 0
      flash.data$Z              <- 1L * is.na(data)
    } else {
      flash.data$Z <- 1
    }

    if (!is.null(S)) {
      if (length(S) == 1) {
        tau.dim <- 0
      } else {
        dims.must.match(data, S)
      }
      tau <- 1 / S^2
    } else if (!is.null(tau)) {
      dims.must.match(data, tau, tau.dim)
    }
    flash.data$tau     <- tau
    flash.data$tau.dim <- tau.dim
  }

  return(flash.data)
}

must.not.supply.both <- function(arg1, arg2) {
  if (!is.null(arg1) && !is.null(arg2))
    stop()
}

must.supply.neither.or.both <- function(arg1, arg2) {
  if ((is.null(arg1) && !is.null(arg2))
      || (!is.null(arg1) && is.null(arg2)))
    stop()
}

must.be.supported.data.type <- function(X,
                                       allow.null = TRUE,
                                       allow.scalar = FALSE,
                                       allow.vector = FALSE) {
  if (!(is.matrix(X)
        || is(X, "Matrix")
        || (is.array(X) && length(dim(X)) == 3)
        || (allow.null && is.null(X))
        || (allow.scalar && is.vector(X) && length(X) == 1)
        || (allow.vector && is.vector(X))))
    stop()
}

must.be.valid.integer <- function(x, lower = NULL, upper = NULL, allow.null = TRUE) {
  if (is.null(x)) {
    if (!allow.null)
      stop(invalid.arg.error(x))
  } else if (!(is.numeric(x)
               && is.finite(x)
               && as.integer(x) == x
               && (is.null(lower) || x >= lower)
               && (is.null(upper) || x <= upper)))
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
