#' Initialize a flash factor
#'
#' Initializes a new flash factor. Custom initialization functions can also be
#'   used. They should accept a single parameter \code{flash} and
#'   output a list consisting of two vectors (which will be interpreted as a
#'   rank-one matrix).
#'
#' @param flash A \code{flash.fit} object.
#'
#' @param tol Convergence tolerance.
#'
#' @param maxiter Maximum number of iterations.
#'
#' @param seed Since initialization is random, a default seed is set for
#'   reproducibility.
#'
#' @seealso init.fn.softImpute
#'
#' @examples
#'
#' # Change the default initialization maxiter
#' my.init.fn <- function(flash) {
#'   return(init.fn.default(flash, maxiter = 500))
#' }
#'
#' fl <- flash.init(gtex) %>% flash.add.greedy(init.fn = my.init.fn)
#'
#' # A wrapper to function nnmf in package NNLM
#' nnmf.init.fn <- function(flash) {
#'   R <- flash$Y - flashier:::lowrank.expand(flash$EF)
#'   res <- NNLM::nnmf(R, verbose = FALSE)
#'   return(list(as.vector(res$W), as.vector(res$H)))
#' }
#'
#' fl.nnmf <- flash.init(gtex) %>%
#'   flash.add.greedy(ebnm.fn = ebnm::ebnm_unimodal_nonnegative,
#'                    init.fn = nnmf.init.fn)
#'
#' @export
#'
init.fn.default <- function(flash,
                            tol = 1 / max(get.dims(flash)),
                            maxiter = 100,
                            seed = 666) {
  return(
    init.fn.constrained(
      flash,
      mode.signs = rep(0, get.dim(flash)),
      tol = tol,
      maxiter = maxiter,
      seed = seed
    )
  )
}

init.fn.constrained <- function(flash,
                              mode.signs,
                              tol = 1 / max(get.dims(flash)),
                              maxiter = 100,
                              seed = 666) {
  set.seed(seed)
  EF <- r1.random(get.dims(flash), mode.signs)

  update.order <- 1:get.dim(flash)
  # Nonnegative/nonpositive dimensions are updated last.
  signed.dims <- which(mode.signs %in% c(-1, 1))
  if (length(signed.dims) > 0) {
    which.signed <- which(update.order %in% signed.dims)
    update.order <- c(update.order[-which.signed], update.order[which.signed])
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.EF <- EF
    EF <- update.init.EF(EF, flash, update.order, mode.signs)
    max.chg <- calc.max.abs.chg(EF, old.EF)
  }

  # Scale EF so that values aren't too different one dimension from another.
  EF <- scale.EF(EF)

  return(EF)
}

update.init.EF <- function(EF, flash, update.order, mode.signs) {
  if (is.null(mode.signs))
    mode.signs <- rep(0, get.dim(flash))

  for (n in update.order) {
    sign <- mode.signs[n]
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

#' Initialize a flash factor using softImpute
#'
#' Initializes a new flash factor using \code{softImpute::softImpute}. This
#'   is slower for very large data matrices, but often yields better results
#'   when there is missing data.
#'
#' @inheritParams init.fn.default
#'
#' @export
#'
init.fn.softImpute <- function(flash,
                               tol = 1 / max(get.dims(flash)),
                               maxiter = 100,
                               seed = 666) {
  set.seed(seed)

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
