update.existing.factor <- function(flash, k, iter, verbose.lvl) {
  old.factor <- extract.factor(flash, k)

  if (!is.zero(old.factor)) {
    factor <- update.factor(old.factor, flash)

    if (get.obj(factor) > get.obj(flash) || !is.obj.valid(flash, factor)) {
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
  factor$obj     <- get.obj(flash)
  factor$is.zero <- is.zero(flash, k)

  factor$fix.dim <- get.fix.dim(flash, k)
  if (!is.null(factor$fix.dim))
    factor$idx.subset <- get.unfixed.idx(flash, k)

  return(factor)
}

alter.existing.factor <- function(flash, factor) {
  k <- get.k(factor)

  flash <- update.residuals(flash, factor)
  flash <- set.EFk(flash, k, get.EF(factor))
  flash <- set.EF2k(flash, k, get.EF2(factor))
  flash <- set.KLk(flash, k, get.KL(factor))
  flash <- set.gk(flash, k, get.g(factor))
  flash <- set.tau(flash, get.tau(factor))
  flash <- set.obj(flash, get.obj(factor))
  if (is.tau.lowrank(flash)) {
    flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
    flash <- set.est.tau(flash, get.est.tau(factor))
  }

  if (is.zero(factor))
    flash <- set.to.zero(flash, k)
  if (is.valid(factor))
    flash <- set.to.valid(flash, k)

  return(flash)
}
