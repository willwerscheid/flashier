set.flash.data <- function(data, S = NULL, S.dim = NULL, var.type = NULL) {
  if (is(data, "flash.data")) {
    if (is.null(data$given.tau)) {
      S <- sqrt(data$given.S2)
    } else {
      S <- 1 / sqrt(data$given.tau)
    }
    S.dim <- data$given.tau.dim

    if ((use.S2(S, S.dim, var.type) && is.null(data$S2.given))
        || (!use.S2(S, S.dim, var.type) && !is.null(data$S2.given))) {
      warning("Data has not been set correctly for the requested variance ",
              "type. Resetting.")
      data <- data$Y
    } else {
      return(data)
    }
  }

  must.be.supported.data.type(data, allow.null = FALSE)
  must.be.supported.data.type(S, allow.vector = TRUE)
  must.be.integer(S.dim, lower = 0, upper = length(dim(data)))
  must.be.valid.var.type(var.type, length(dim(data)))

  flash.data <- list()
  flash.data$Y <- data
  any.missing <- anyNA(data)
  if (any.missing) {
    flash.data$Y[is.na(data)] <- 0
    flash.data$Z <- 1L * !is.na(data)
  } else {
    flash.data$Z <- 1
  }

  # Set S.dim.
  if (is.vector(S) && is.null(S.dim)) {
    if (length(S) == 1) {
      S.dim <- 0
    } else {
      S.dim <- which(length(S) == dim(data))
      if (length(S.dim) == 0)
        stop("S was interpreted as a vector, but couldn't be aligned ",
             "with the data.")
      if (length(S.dim) > 1)
        stop("S could not be unambiguously interpreted. Set data using ",
             "set.flash.data with S.dim specified.")
    }
  } else {
    dims.must.match(data, S, S.dim)
  }

  # S2 is stored for "noisy" variance structures. Otherwise tau is stored.
  if (use.S2(S, S.dim, var.type)) {
    S2 <- S^2
    if (is.vector(S2)) {
      # Convert S2 to a matrix or array.
      each <- prod(c(1, dim(data)[1:length(dim(data)) < S.dim]))
      S2 <- array(rep(S2, each = each), dim = dim(data))
    }
    S2[is.na(data)] <- Inf
    flash.data$given.S2 <- S2
  } else if (!is.null(S)) {
    tau <- 1 / S^2
    if (is.null(S.dim))
      tau[is.na(data)] <- 0
    flash.data$given.tau <- tau
    flash.data$given.tau.dim <- S.dim
  }

  class(flash.data) <- "flash.data"

  return(flash.data)
}

use.S2 <- function(S, S.dim, var.type) {
  return(!is.null(S)
         && (length(var.type) == 1)
         && (is.null(S.dim) || (S.dim > 0 && !(S.dim == var.type))))
}

