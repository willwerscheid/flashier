# Set up the data for flash. A "flash.data" object always contains fields "Y"
#   (the data, with NAs replaced by zeros), "Z" (a matrix or array that has
#   zeros where data is missing and ones where data is nonmissing; if no data
#   is missing, then Z = 1), and "est.tau.dim" (the variance type). If standard
#   errors are provided, then (depending on the variance type) the object will
#   also contain fields "given.S2" or "given.tau" and "given.tau.dim".
#
set.flash.data <- function(data, S = NULL, S.dim = NULL, var.type = NULL) {
  flash.data <- list(t.init = Sys.time(), t.final = Sys.time())

  data <- handle.data(data)

  must.be.supported.S.type(S)
  must.be.compatible.data.types(data, S)

  data.dim <- get.data.dims(data)
  must.be.integer(S.dim, lower = 0, upper = length(data.dim), allow.null = TRUE)
  must.be.valid.var.type(var.type, length(data.dim))

  # Set Y, Z, and est.tau.dim.
  flash.data$Y <- data
  if (anyNA(data)) {
    flash.data$Y[is.na(data)] <- 0
    flash.data$Z <- 1L * !is.na(data)
  } else {
    flash.data$Z <- 1
  }
  flash.data$est.tau.dim <- var.type

  must.not.have.zero.slices(flash.data$Y, S, var.type)

  # Set S.dim.
  if (is.null(S.dim)) {
    S.dim <- infer.S.dim(S, data.dim)
  } else {
    dims.must.match(data, S, S.dim)
  }

  if (is.S2.stored(S, S.dim, var.type)) {
    # TODO: estimate.noisy.kron.tau cannot yet handle missing data.
    if (anyNA(data) && (length(var.type) > 1)) {
      stop("The noisy Kronecker variance structure has not yet been implemented ",
           "for missing data.")
    }

    S2 <- S^2

    if (is.vector(S)) {
      # TODO: estimate.noisy.tau and estimate.noisy.kron.tau cannot yet handle
      #   vector-valued S.
      each <- prod(c(1, data.dim[1:length(data.dim) < S.dim]))
      S2 <- array(rep(S2, each = each), dim = data.dim)
    }

    S2[is.na(data)] <- Inf

    flash.data$given.S2 <- S2
  } else if (!is.null(S)) {
    tau <- 1 / S^2

    if (!is.vector(S)) {
      tau[is.na(data)] <- 0
    }

    flash.data$given.tau <- tau
    flash.data$given.tau.dim <- S.dim
  }

  class(flash.data) <- c("flash.data", "list")

  return(flash.data)
}

handle.data <- function(X) {
  # Low-rank and low-rank plus sparse representations as named lists:
  if (is.list(X) &&
      any(c("u", "U") %in% names(X)) &&
      any(c("v", "V") %in% names(X))) {
    X <- handle.udvs(X)
  }

  return(X)
}

# Low-rank or low-rank plus sparse representations, of form:
#   X = u diag(d) v' + s.
#
handle.udvs <- function(X) {
  names(X) <- tolower(names(X))

  # Fields u and v (or U and V) are required. d and s (or D and S) are optional.
  u <- X$u
  v <- X$v
  d <- X$d
  s <- X$s

  # Check u and v.
  if (!(is.matrix(u) || inherits(u, "Matrix"))
      || !(is.matrix(v) || inherits(v, "Matrix"))) {
    stop("Data fields u and v must be matrices (or sparse matrices of class ",
         "\"Matrix\").")
  }
  if (anyNA(u) || anyNA(v)) {
    stop("Data matrices u and v must not have missing data.")
  }

  if (!identical(ncol(u), ncol(v))) {
    # Silently handle the case where v^T is provided (rather than v).
    if (identical(ncol(u), nrow(v))) {
      v <- t(v)
    } else {
      stop("Data matrices u and v are of incompatible dimensions.")
    }
  }

  K <- ncol(u)

  # Check d.
  if (!is.null(d)) {
    if (!is.numeric(d) || (length(d) < K)) {
      stop("If data field d is provided, it must be a vector of length K.")
    }
    d <- d[1:K]
    u <- t(t(u) * sqrt(d))
    v <- t(t(v) * sqrt(d))
  }

  udv <- list(u, v)
  class(udv) <- c("lowrank", "list")

  # Check s.
  if (is.null(s)) {
    ret <- udv
  } else {
    if (!(is.matrix(s) || inherits(s, "Matrix"))) {
      stop("If data field s is provided, it must be a matrix or object of ",
           "class \"Matrix\".")
    }
    if (!identical(nrow(u), nrow(s)) || !identical(nrow(v), ncol(s))) {
      stop("The dimensions of data matrix s do not match the dimensions of ",
           "uv^T.")
    }
    ret <- list("lowrank" = udv, "sparse" = s)
    class(ret) <- c("lrps", "list")
  }

  return(ret)
}

# If S is a vector, then attempt to infer its mode (e.g., row-wise, column-wise).
#
infer.S.dim <- function(S, data.dim) {
  if (!is.vector(S)) {
    S.dim <- NULL
  } else if (length(S) == 1) {
    S.dim <- 0
  } else {
    S.dim <- which(length(S) == data.dim)
    if (length(S.dim) == 0) {
      stop("S was interpreted as a vector, but couldn't be aligned ",
           "with the data.")
    } else if (length(S.dim) > 1) {
      stop("S could not be unambiguously interpreted. Please specify S.dim.")
    }
  }

  return(S.dim)
}

# When given and estimated variances vary along one and the same mode (for
#   example, when row-wise standard errors are provided and additional
#   row-wise variances are to be estimated), then providing standard errors
#   is equivalent to putting lower limits on the estimated standard errors.
#   Otherwise, estimation is a difficult optimization problem. In the former
#   case, it is more convenient to store tau; in the latter, it is better to
#   store S^2.
#
is.S2.stored <- function(S, S.dim, var.type) {
  return(!is.null(S)
         && (length(var.type) > 0)
         && (length(var.type) > 1
             || is.null(S.dim)
             || (S.dim > 0 && !(S.dim == var.type))))
}
