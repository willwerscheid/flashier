update.existing.factor <- function(flash, k, iter, verbose.lvl) {
  old.factor <- extract.factor(flash, k)

  if (!is.zero(old.factor)) {
    factor <- update.factor(old.factor, flash)

    if (get.obj(factor) > get.obj(flash)
        || !is.obj.valid(flash, factor)
        || !identical(get.exclusions(old.factor), get.exclusions(factor))) {
      flash <- alter.existing.factor(flash, factor)
    } else {
      obj.diff <- get.obj(factor) - get.obj(flash)
      report.backfit.obj.decrease(verbose.lvl, obj.diff, k)
    }
  }

  return(flash)
}

extract.factor <- function(flash, k) {
  factor         <- list()
  factor$k       <- k
  factor$EF      <- get.EF.k(flash, k)
  factor$EF2     <- get.EF2.k(flash, k)
  factor$KL      <- get.KL.k(flash, k)
  factor$g       <- get.g.k(flash, k)
  factor$obj     <- get.obj(flash)
  factor$is.zero <- is.zero(flash, k)

  factor$fix.dim <- get.fix.dim(flash, k)
  if (!is.null(factor$fix.dim))
    factor$idx.subset <- get.unfixed.idx(flash, k)

  factor$exclusions <- get.exclusions(flash, k)

  return(factor)
}

alter.existing.factor <- function(flash, factor) {
  k <- get.k(factor)

  if (uses.R(flash))
    flash <- update.R(flash, factor)

  flash <- set.EFk(flash, k, get.EF(factor))
  flash <- set.EF2k(flash, k, get.EF2(factor))
  flash <- set.KLk(flash, k, get.KL(factor))
  flash <- set.gk(flash, k, get.g(factor))
  flash <- set.tau(flash, get.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  if (is.tau.simple(flash)) {
    flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
    flash <- set.est.tau(flash, get.est.tau(factor))
  }

  if (is.zero(factor))
    flash <- set.to.zero(flash, k)
  if (is.valid(factor))
    flash <- set.to.valid(flash, k)
  flash <- add.exclusions(flash, get.exclusions(factor), k)

  return(flash)
}

update.factors.parallel <- function(flash, kset) {
  factors <- lapply(1:get.n.factors(flash), function(k) {
    factor <- extract.factor(flash, k)
    if (k %in% kset) {
      factor <- update.factor(factor, flash, update.tau = FALSE)
      factor <- set.to.valid(factor)
    }
    return(factor)
  })

  flash <- alter.factors.parallel(flash, factors)

  return(flash)
}

alter.factors.parallel <- function(flash, factors) {
  d <- get.dim(flash)
  EF  <- lapply(1:d, function(n) {sapply(factors, function(f) get.EF(f)[[n]])})
  EF2 <- lapply(1:d, function(n) {sapply(factors, function(f) get.EF2(f)[[n]])})
  class(EF)  <- "lowrank"
  class(EF2) <- "lowrank"

  flash <- set.EF(flash, EF)
  flash <- set.EF2(flash, EF2)
  flash <- set.KL(flash, lapply(1:d, function(n) {
    sapply(factors, function(f) get.KL(f)[[n]])
  }))
  flash <- set.g(flash, lapply(factors, getElement, "g"))
  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))
  flash <- set.is.valid(flash, sapply(factors, getElement, "is.valid"))
  flash <- set.is.zero(flash, sapply(factors, getElement, "is.zero"))

  return(flash)
}
