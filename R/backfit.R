# TODO: delete backfit and backfit.once

backfit <- function(flash, kset, shuffle.kset = FALSE, tol = 1e-2) {
  flash <- optimize.it(flash,
                       update.fn = backfit.once,
                       update.args = list(kset = kset,
                                          shuffle.kset = shuffle.kset),
                       obj.fn = calc.obj.diff,
                       tol = tol)
  return(flash)
}

backfit.once <- function(flash, kset, shuffle.kset = FALSE) {
  if (shuffle.kset)
    kset <- sample(kset)

  for (k in kset)
    flash <- update.kth.factor(flash, k)

  return(flash)
}

update.kth.factor <- function(flash, k) {
  factor <- extract.factor(flash, k)

  if (!is.zero(factor)) {
    factor <- update.factor(factor, flash)
    if (get.obj(factor) > get.obj(flash))
      flash <- alter.existing.factor(flash, factor)
    # TODO: what if objective decreases?
  }

  return(flash)
}

extract.factor <- function(flash, k) {
  factor         <- list()
  factor$EF      <- get.EFk(flash, k)
  factor$EF2     <- get.EF2k(flash, k)
  factor$est.tau <- get.est.tau(flash)
  factor$is.zero <- is.zero(flash, k)
  factor$k       <- k

  return(factor)
}

alter.existing.factor <- function(flash, factor) {
  k <- get.k(factor)

  if (uses.R(flash)) {
    new.EF <- as.lowrank(get.EF(factor))
    old.EF <- as.lowrank(get.EFk(flash, k))
    EF.delta.mat <- lowrank.delta.mat(new.EF, old.EF)
    flash$R <- flash$R - get.nonmissing(flash) * lowrank.expand(EF.delta.mat)
  }

  flash <- set.EFk(flash, k, get.EF(factor))
  flash <- set.EF2k(flash, k, get.EF2(factor))
  flash <- set.KLk(flash, k, get.KL(factor))
  flash <- set.gk(flash, k, get.g(factor))
  flash <- set.R2(flash, get.R2(flash) + get.delta.R2(factor))
  flash <- set.est.tau(flash, get.est.tau(factor))
  flash <- set.obj(flash, get.obj(factor))

  if (is.zero(factor))
    flash <- set.to.zero(flash, k)

  return(flash)
}
