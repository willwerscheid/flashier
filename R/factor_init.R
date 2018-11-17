init.factor <- function(flash) {
  factor          <- list()
  factor$EF       <- init.r1(flash)
  factor$EF2      <- r1.square(factor$EF)
  factor$delta.R2 <- calc.delta.R2(factor, flash)
  factor$est.tau  <- calc.est.tau(flash, factor$delta.R2)
  factor$obj      <- -Inf # not yet a valid factor
  factor$is.zero  <- FALSE

  return(factor)
}

init.r1 <- function(flash, tol = 1e-3) {
  R  <- get.R(flash)
  Y  <- get.Y(flash)
  Z  <- get.nonmissing(flash)
  EF <- get.EF(flash)

  fix.dim     <- get.next.fix.dim(flash)
  fix.idx     <- get.next.fix.idx(flash)
  fix.vals    <- get.next.fix.vals(flash)
  nonneg.dims <- get.next.nonneg.dims(flash)

  r1 <- r1.random(get.dims(flash), nonneg.dims)
  if (!is.null(fix.dim))
    r1[[fix.dim]][fix.idx] <- fix.vals

  update.order <- 1:length(r1)
  # Nonnegative dimensions are updated last so that they stay nonnegative:
  if (!is.null(nonneg.dims))
    update.order <- c(update.order[-nonneg.dims], update.order[nonneg.dims])

  subset.data <- list()
  if (!is.null(fix.dim)) {
    subset.data$idx.subset <- setdiff(1:length(r1[[fix.dim]]), fix.idx)

    # If a dimension is entirely fixed, then it never needs to be updated:
    if (length(subset.data$idx.subset) == 0) {
      update.order <- setdiff(update.order, fix.dim)
    # Otherwise, updates are performed using subsetted matrices:
    } else {
      subset.data$R.subset  <- fullrank.subset(R, fix.dim, subset.data$idx.subset)
      subset.data$Y.subset  <- fullrank.subset(Y, fix.dim, subset.data$idx.subset)
      subset.data$Z.subset  <- fullrank.subset(Z, fix.dim, subset.data$idx.subset)
      subset.data$EF.subset <- lowrank.subset(EF, fix.dim, subset.data$idx.subset)
    }
  }

  r1 <- optimize.it(r1,
                    update.fn = update.r1,
                    update.args = list(flash = flash,
                                       update.order = update.order,
                                       fix.dim = fix.dim,
                                       nonneg.dims = nonneg.dims,
                                       subset.data = subset.data),
                    obj.fn = calc.max.chg.r1,
                    tol = tol)

  return(r1)
}

update.r1 <- function(r1, flash, update.order, fix.dim, nonneg.dims, subset.data) {
  for (n in update.order) {
    is.fixed  <- (identical(n, fix.dim))
    is.nonneg <- (n %in% nonneg.dims)

    r1 <- update.r1.one.n(r1, n, flash, is.fixed, is.nonneg, subset.data)
  }

  return(r1)
}

update.r1.one.n <- function(r1, n, flash, is.fixed, is.nonneg, subset.data) {
  if (is.fixed) {
    R  <- subset.data$R.subset
    Y  <- subset.data$Y.subset
    Z  <- subset.data$Z.subset
    EF <- subset.data$EF.subset
  } else {
    R  <- get.R(flash)
    Y  <- get.Y(flash)
    Z  <- get.nonmissing(flash)
    EF <- get.EF(flash)
  }

  if (uses.R(flash)) {
    new.vals <- (nmode.prod.r1(R, r1[-n], n)
                 / nmode.prod.r1(Z, r1.square(r1[-n]), n))
  } else {
    new.vals <- ((nmode.prod.r1(Y, r1[-n], n)
                  - premult.nmode.prod.r1(Z, EF, r1[-n], n))
                 / nmode.prod.r1(Z, r1.square(r1[-n]), n))
  }

  if (is.nonneg)
    new.vals <- pmax(new.vals, 0)

  if (is.fixed) {
    r1[[n]][subset.data$idx.subset] <- new.vals
  } else {
    r1[[n]] <- new.vals
  }

  return(r1)
}
