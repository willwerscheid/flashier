#   A custom initialization function may also be used. It should accept
# parameters flash, tol, and maxiter and output a list of vectors (which will
# be interpreted as an "r1" object). For example, a wrapper to function nnmf
# in package NNLM can be written as follows:
#
# nnmf.init.fn <- function(flash, tol, maxiter) {
#   res <- NNLM::nnmf(flash$Y,
#                     init = list(W0 = flash$EF[[1]],
#                                 H0 = t(flash$EF[[2]])),
#                     rel.tol = tol,
#                     max.iter = maxiter,
#                     verbose = FALSE)
#   return(list(as.vector(res$W), as.vector(res$H)))
# }

init.factor <- function(flash, init.fn, tol, maxiter) {
  factor <- list()
  factor$is.fixed <- is.next.fixed(flash)

  if (is.null(init.fn)) {
    factor$EF <- init.next.EF(flash, tol, maxiter)
  } else {
    factor$EF <- do.call(init.fn, list(flash, tol, maxiter))
    class(factor$EF) <- "r1"
  }

  factor$EF2      <- r1.square(factor$EF)
  factor$KL       <- rep(0, get.dim(flash))
  factor          <- update.tau(factor, flash)
  factor$obj      <- calc.obj(flash, factor)
  factor$is.valid <- FALSE
  factor$is.zero  <- FALSE

  return(factor)
}

init.next.EF <- function(flash, tol, maxiter) {
  next.k    <- get.next.k(flash)
  fix.dim   <- get.fix.dim(flash, next.k)
  fix.idx   <- get.fix.idx(flash, next.k)
  fix.vals  <- get.fix.vals(flash, next.k)
  dim.signs <- get.dim.signs(flash, next.k)

  EF <- r1.random(get.dims(flash), dim.signs)
  if (!is.null(fix.dim)) {
    EF[[fix.dim]][fix.idx] <- fix.vals
    # The non-fixed values should be similar in magnitude to the fixed ones.
    mean.fix2 <- mean(fix.vals^2)
    if (mean.fix2 != 0)
      EF[[fix.dim]][-fix.idx] <- sqrt(mean.fix2) * EF[[fix.dim]][-fix.idx]
  }

  update.order <- 1:get.dim(flash)
  # Nonnegative/nonpositive dimensions are updated later.
  signed.dims <- which(dim.signs %in% c(-1, 1))
  if (length(signed.dims) > 0) {
    which.signed <- which(update.order %in% signed.dims)
    update.order <- c(update.order[-which.signed], update.order[which.signed])
  }
  # And fixed dimensions are updated last.
  if (!is.null(fix.dim)) {
    which.fixed <- which(update.order == fix.dim)
    update.order <- c(update.order[-which.fixed], which.fixed)
  }

  subset.data <- NULL
  if (!is.null(fix.dim)) {
    idx.subset <- get.next.unfixed.idx(flash)
    # If a dimension is entirely fixed, then it never needs to be updated.
    if (length(idx.subset) == 0) {
      update.order <- setdiff(update.order, fix.dim)
      # Otherwise, updates are performed using subsetted matrices.
    } else {
      subset.data <- get.subset.data(flash, fix.dim, idx.subset)
    }
  }

  max.chg <- Inf
  iter <- 0
  while (max.chg > tol && iter < maxiter) {
    iter <- iter + 1
    old.EF <- EF
    EF <- update.init.EF(EF, flash, update.order, fix.dim, dim.signs,
                         subset.data)
    max.chg <- calc.max.abs.chg(EF, old.EF)
  }

  # Scale EF so that values aren't too different one dimension from another.
  if (is.null(fix.dim))
    EF <- scale.EF(EF)

  return(EF)
}

update.init.EF <- function(EF, flash, update.order, fix.dim, dim.signs,
                           subset.data) {
  if (is.null(dim.signs))
    dim.signs <- rep(0, get.dim(flash))

  for (n in update.order) {
    is.fixed <- n %in% fix.dim
    sign <- dim.signs[n]
    EF <- update.init.EF.one.n(EF, n, flash, is.fixed, sign, subset.data)
  }

  return(EF)
}

update.init.EF.one.n <- function(EF, n, flash, is.fixed, sign, subset.data) {
  if (is.fixed) {
    R        <- subset.data$R.subset
    Y        <- subset.data$Y.subset
    Z        <- subset.data$Z.subset
    flash.EF <- subset.data$EF.subset
  } else {
    R        <- get.R(flash)
    Y        <- get.Y(flash)
    Z        <- get.nonmissing(flash)
    flash.EF <- get.EF(flash)
  }

  if (uses.R(flash)) {
    new.vals <- (nmode.prod.r1(R, EF[-n], n)
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  } else {
    new.vals <- ((nmode.prod.r1(Y, EF[-n], n)
                  - premult.nmode.prod.r1(Z, flash.EF, EF[-n], n))
                 / nmode.prod.r1(Z, r1.square(EF[-n]), n))
  }

  if (sign == 1)
    new.vals <- pmax(new.vals, 0)
  if (sign == -1)
    new.vals <- pmin(new.vals, 0)

  if (is.fixed) {
    EF[[n]][subset.data$idx.subset] <- new.vals
  } else {
    EF[[n]] <- new.vals
  }

  return(EF)
}

scale.EF <- function(EF) {
  norms <- lapply(EF, function(x) {sqrt(sum(x^2))})
  EF <- mapply(`/`, EF, norms, SIMPLIFY = FALSE)
  EF <- lapply(EF, `*`, prod(unlist(norms))^(1/length(EF)))
  class(EF) <- "r1"
  return(EF)
}
