#' Set data for flash
#'
#' Converts matrices or arrays of observations and standard errors into the
#' form that will be used by \code{\link{flashier}}. When it is necessary to be
#' parsimonious with memory, one can call \code{set.flash.data} and then remove
#' the original objects from memory. Otherwise, one should simply pass the
#' original objects to \code{flashier} as is.
#'
#' @inheritParams flashier
#'
#' @param S.dim The dimension along which \code{S} lies when \code{S} is a
#'   vector. Only necessary when it cannot be inferred from the data (when,
#'   for example, \code{data} is a square matrix).
#'
#' @export

set.flash.data <- function(data, S = NULL, S.dim = NULL, var.type = NULL) {
  # If data is a flash.data object, check that it has been set correctly.
  if (inherits(data, "flash.data")) {
    if (!is.null(data$given.S2)) {
      S <- data$given.S2
    } else if (!is.null(data$given.tau)) {
      S <- data$given.tau
    }
    S.dim <- data$given.tau.dim

    if ((use.S2(S, S.dim, var.type) && is.null(data$given.S2))
        || (!use.S2(S, S.dim, var.type) && !is.null(data$given.S2))) {
      warning("Data has not been set correctly for the requested variance ",
              "type. Resetting.")
      data <- data$Y
      if (!is.null(data$given.S2)) {
        S <- sqrt(S)
      } else if (!is.null(data$given.tau)) {
        S <- 1 / sqrt(S)
      }
    } else {
      return(data)
    }
  }

  # If data is a list, attempt to interpret it as a low-rank representation
  #   of the data. Must have fields d, u, and v.
  if (is.list(data)) {
    error.msg <- paste("Data is a list but could not be interpreted as a",
                       "low-rank representation.")

    if (is.null(data$d) || is.null(data$u) || is.null(data$v)) {
      stop(error.msg)
    }

    if (!is.matrix(data$u) || !is.matrix(data$v)
        || !identical(ncol(data$u), ncol(data$v))) {
      stop(error.msg)
    }

    K <- ncol(data$u)

    if (!is.numeric(data$d) || (length(data$d) < K)) {
      stop(error.msg)
    }

    LR <- list(data$u, data$v)
    D  <- data$d[1:K]

    LR[[1]] <- t(t(LR[[1]]) * sqrt(D))
    LR[[2]] <- t(t(LR[[2]]) * sqrt(D))

    if (any(sapply(LR, anyNA)))
      stop("If a low-rank representation of the data is used, then no data can",
           " be missing.")

    class(LR) <- "lowrank"
    data <- LR
    dim.data <- sapply(LR, ncol)
  } else {
    dim.data <- dim(data)
  }

  must.be.supported.data.type(data, allow.null = FALSE, allow.lowrank = TRUE)
  must.be.supported.data.type(S, allow.vector = TRUE)
  must.be.compatible.data.types(data, S)
  must.be.integer(S.dim, lower = 0, upper = length(dim.data))
  must.be.valid.var.type(var.type, length(dim.data))

  # Set Y and Z.
  flash.data <- list()
  flash.data$Y <- data
  any.missing <- anyNA(data)
  if (any.missing) {
    flash.data$Y[is.na(data)] <- 0
    flash.data$Z <- 1L * !is.na(data)
  } else {
    flash.data$Z <- 1
  }
  must.not.have.zero.slices(flash.data$Y)

  # Set S.dim.
  if (is.vector(S) && is.null(S.dim)) {
    if (length(S) == 1) {
      S.dim <- 0
    } else {
      S.dim <- which(length(S) == dim.data)
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
      each <- prod(c(1, dim.data[1:length(dim.data) < S.dim]))
      S2 <- array(rep(S2, each = each), dim = dim.data)
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
         && (length(var.type) > 0)
         && (length(var.type) > 1
             || is.null(S.dim)
             || (S.dim > 0 && !(S.dim == var.type))))
}
