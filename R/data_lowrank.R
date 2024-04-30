# Low-rank matrix representations: X = UV' or X_ijm = \sum_k U_ik V_jk W_mk.
#
# Internally, an object of class "lowrank" is an unnamed list consisting of
# two matrices (for matrix data) or three matrices (for tensor data). All
# matrices should have the same number of columns; e.g., if X is n x p, the
# lowrank object would consist of a n x k matrix and a p x k matrix.
#
#' @export
as.fullrank.lowrank <- function(x, ...) {
  return(lowrank.expand(x))
}

#' @export
get.data.dims.lowrank <- function(x, ...) {
  return(sapply(x, nrow))
}

#' @export
get.data.dimnames.lowrank <- function(x, ...) {
  return(lapply(x, rownames))
}

#' @export
nmode.prod.vec.lowrank <- function(x, v, n, ...) {
  if (length(x) == 2) {
    if (identical(v, 1))
      u <- colSums(x[[n]])
    else
      u <- t(v %*% x[[n]])
    return(as.vector(x[[-n]] %*% u))
  } else if (length(x) == 3) {
    if (identical(v, 1))
      u <- colSums(x[[n]])
    else
      u <- as.vector(v %*% x[[n]])

    retval <- x[-n]
    smaller.dim <- which.min(c(nrow(retval[[1]]), nrow(retval[[2]])))
    retval[[smaller.dim]] <- t(t(retval[[smaller.dim]]) * u)
    class(retval) <- "lowrank"
    return(retval)
  } else {
    stop(nmode.ops.error)
  }
}

#' @export
premult.nmode.prod.r1.lowrank <- function(x, lr, r1, n, ...) {
  if (is.null(lr))
    return(0)

  if (inherits(lr, "lowrank")) {
    if (length(lr) == 2) {
      u <- unlist(r1) * lr[[-n]]
      if (n == 1)
        return(rowSums(lr[[1]] * (x[[1]] %*% crossprod(x[[2]], u))))
      if (n == 2)
        return(rowSums(lr[[2]] * (x[[2]] %*% t(crossprod(u, x[[1]])))))
    } else if (length(lr) == 3) {
      return(premult.nmode.prod.r1.default(x, lr, r1, n, ...))
    } else {
      stop(nmode.ops.error)
    }
  }

  # else lr is full-rank:
  return(premult.nmode.prod.r1.lowrank(lr, x, r1, n, ...))
}

#' @export
sq.nmode.prod.r1.lowrank <- function(x, r1, n, ...) {
  if (length(x) == 2 && ncol(x[[1]])^2 > prod(get.data.dims(x)))
    return(nmode.prod.r1(lowrank.expand(x)^2, r1, n))
  return(premult.nmode.prod.r1(x, x, r1, n))
}

# Not currently used -----

# Only needed for parallel backfits:
#' @export
nmode.prod.lowrank.lowrank <- function(x, Y, n, ...) {
  if (length(x) == 2) {
    if (n == 1)
      return(x[[1]] %*% (t(x[[2]]) %*% Y[[1]]))
    if (n == 2)
      return(x[[2]] %*% (t(x[[1]]) %*% Y[[1]]))
  } else if (length(x) == 3) {
    stop("N-mode products of tensors with matrices have not yet been ",
         "implemented.")
  } else {
    stop(nmode.ops.error)
  }
}
