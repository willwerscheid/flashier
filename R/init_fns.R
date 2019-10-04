#   Custom initialization functions may also be used. They should accept
# parameters flash and dim.signs and output a list of vectors (which will
# be interpreted as an "r1" object). For example, a wrapper to function nnmf
# in package NNLM can be written as follows:
#
# nnmf.init.fn <- function(flash, dim.signs) {
#   R <- flash$Y - flashier:::lowrank.expand(flash$EF)
#   res <- NNLM::nnmf(R, rel.tol = tol, max.iter = maxiter, verbose = FALSE)
#   return(list(as.vector(res$W), as.vector(res$H)))
# }

# TODO document

#' @export
#'
init.fn.softImpute <- function(flash,
                               dim.signs,
                               tol = 1 / max(get.dims(flash)),
                               maxiter = 100) {
  if (get.dim(flash) > 2)
    stop("softImpute cannot be used with tensors.")

  if (inherits(get.Y(flash), "lowrank"))
    stop("softImpute cannot be used with low-rank matrix representations.")

  if (uses.R(flash)) {
    R <- get.R(flash)
  } else {
    R <- get.Y(flash) - lowrank.expand(get.EF(flash))
  }

  if (any.missing(flash)) {
    R[get.nonmissing(flash) == 0] <- NA
  }

  suppressWarnings({
    si.res <- softImpute::softImpute(R, rank.max = 1, type = "als", lambda = 0)
  })

  EF <- list(si.res$u * sqrt(si.res$d), si.res$v * sqrt(si.res$d))
  class(EF) <- "r1"

  return(EF)
}

#' @export
#'
init.fn.default <- function(flash,
                            dim.signs,
                            tol = 1 / max(get.dims(flash)),
                            maxiter = 100) {
  EF <- r1.random(get.dims(flash), dim.signs)

  update.order <- 1:get.dim(flash)
  # Nonnegative/nonpositive dimensions are updated last.
  signed.dims <- which(dim.signs %in% c(-1, 1))
  if (length(signed.dims) > 0) {
    which.signed <- which(update.order %in% signed.dims)
    update.order <- c(update.order[-which.signed], update.order[which.signed])
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.EF <- EF
    EF <- update.init.EF(EF, flash, update.order, dim.signs)
    max.chg <- calc.max.abs.chg(EF, old.EF)
  }

  # Scale EF so that values aren't too different one dimension from another.
  EF <- scale.EF(EF)

  return(EF)
}

update.init.EF <- function(EF, flash, update.order, dim.signs) {
  if (is.null(dim.signs))
    dim.signs <- rep(0, get.dim(flash))

  for (n in update.order) {
    sign <- dim.signs[n]
    EF <- update.init.EF.one.n(EF, n, flash, sign)
  }

  return(EF)
}

update.init.EF.one.n <- function(EF, n, flash, sign) {
  R        <- get.R(flash)
  Y        <- get.Y(flash)
  Z        <- get.nonmissing(flash)
  flash.EF <- get.EF(flash)

  if (uses.R(flash)) {
    new.vals <- (nmode.prod.r1(R, EF[-n], n)
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  } else {
    new.vals <- ((nmode.prod.r1(Y, EF[-n], n)
                  - premult.nmode.prod.r1(Z, flash.EF, EF[-n], n))
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  }

  new.vals[is.na(new.vals)] <- 0

  if (sign == 1)
    new.vals <- pmax(new.vals, 0)
  if (sign == -1)
    new.vals <- pmin(new.vals, 0)

  EF[[n]] <- new.vals

  return(EF)
}

scale.EF <- function(EF) {
  norms <- lapply(EF, function(x) {sqrt(sum(x^2))})

  if (all(unlist(norms) > 0)) {
    EF <- mapply(`/`, EF, norms, SIMPLIFY = FALSE)
    EF <- lapply(EF, `*`, prod(unlist(norms))^(1/length(EF)))
  } else {
    warning("Fitting stopped after the initialization function failed to find",
            " a non-zero factor.")
    EF <- lapply(EF, `*`, 0)
  }
  class(EF) <- "r1"

  return(EF)
}
