init.factor <- function(flash, tol = 1e-2) {
  factor             <- list()
  factor$is.fixed    <- is.next.fixed(flash)
  factor$EF          <- init.next.EF(flash, tol)
  factor$EF2         <- r1.square(factor$EF)
  factor$KL          <- rep(0, get.dim(flash))
  factor$delta.R2    <- calc.delta.R2(factor, flash)
  factor$est.tau     <- calc.est.tau(flash, factor$delta.R2)
  factor$obj         <- calc.obj(flash, factor)
  factor$is.valid    <- FALSE
  factor$is.zero     <- FALSE

  return(factor)
}

init.next.EF <- function(flash, tol = 1e-2) {
  fix.dim     <- get.next.fix.dim(flash)
  fix.idx     <- get.next.fix.idx(flash)
  fix.vals    <- get.next.fix.vals(flash)
  nonneg.dims <- get.next.nonneg.dims(flash)

  EF <- r1.random(get.dims(flash), nonneg.dims)
  if (!is.null(fix.dim)) {
    EF[[fix.dim]][fix.idx] <- fix.vals
    # The non-fixed values should be similar in magnitude to the fixed ones:
    if (mean(fix.vals) != 0)
      EF[[fix.dim]][-fix.idx] <- mean(fix.vals) * EF[[fix.dim]][-fix.idx]
  }

  update.order <- 1:get.dim(flash)
  # Nonnegative dimensions are updated later so that they stay nonnegative:
  if (!is.null(nonneg.dims)) {
    which.nonneg <- which(update.order %in% nonneg.dims)
    update.order <- c(update.order[-which.nonneg], update.order[which.nonneg])
  }
  # And fixed dimensions are updated last:
  if (!is.null(fix.dim)) {
    which.fixed  <- which(update.order == fix.dim)
    update.order <- c(update.order[-which.fixed], which.fixed)
  }

  subset.data <- NULL
  if (!is.null(fix.dim)) {
    idx.subset <- get.next.unfixed.idx(flash)
    # If a dimension is entirely fixed, then it never needs to be updated:
    if (length(idx.subset) == 0) {
      update.order <- setdiff(update.order, fix.dim)
    # Otherwise, updates are performed using subsetted matrices:
    } else {
      subset.data <- get.subset.data(flash, fix.dim, idx.subset)
    }
  }

  EF <- optimize.it(EF,
                    update.fn = update.init.EF,
                    update.args = list(flash = flash,
                                       update.order = update.order,
                                       fix.dim = fix.dim,
                                       nonneg.dims = nonneg.dims,
                                       subset.data = subset.data),
                    obj.fn = calc.max.chg.r1,
                    tol = tol)

  return(EF)
}

update.init.EF <- function(EF, flash, update.order, fix.dim, nonneg.dims,
                           subset.data) {
  for (n in update.order) {
    is.fixed  <- n %in% fix.dim
    is.nonneg <- n %in% nonneg.dims
    EF <- update.init.EF.one.n(EF, n, flash, is.fixed, is.nonneg, subset.data)
  }

  return(EF)
}

update.init.EF.one.n <- function(EF, n, flash, is.fixed, is.nonneg,
                                 subset.data) {
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

  if (is.nonneg)
    new.vals <- pmax(new.vals, 0)

  if (is.fixed) {
    EF[[n]][subset.data$idx.subset] <- new.vals
  } else {
    EF[[n]] <- new.vals
  }

  return(EF)
}
