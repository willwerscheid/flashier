# Set up the data for flash. A "flash.data" object always contains fields "Y"
#   (the data, with NAs replaced by zeros), "Z" (a matrix or array that has
#   zeros where data is missing and ones where data is nonmissing; if no data
#   is missing, then Z = 1), and "est.tau.dim" (the variance type). If standard
#   errors are provided, then (depending on the variance type) the object will
#   also contain fields "given.S2" or "given.tau" and "given.tau.dim".
#
set.flash.data <- function(data, S = NULL, S.dim = NULL, var.type = NULL) {
  flash.data <- list(t.init = Sys.time(), t.final = Sys.time())

  must.be.supported.data.type(data, allow.null = FALSE, allow.lowrank = TRUE)
  must.be.supported.data.type(S, allow.vector = TRUE)
  must.be.compatible.data.types(data, S)

  if (is.udv(data)) {
    data <- udv.to.lowrank(data)
  }

  if (inherits(data, "lowrank")) {
    data.dim <- sapply(data, nrow)
  } else {
    data.dim <- dim(data)
  }

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

  must.not.have.zero.slices(flash.data$Y)

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

# Convert a udv-type object to a flashier "lowrank" object.
#
udv.to.lowrank <- function(udv) {
  d <- udv$d[1:ncol(udv$u)]

  LR <- list(t(t(udv$u) * sqrt(d)), t(t(udv$v) * sqrt(d)))
  class(LR) <- c("lowrank", "list")

  return(LR)
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
