nullchk.kth.factor <- function(flash, k) {
  if (!is.valid(flash, k))
    return(flash)

  factor <- zero.factor(flash, k)
  factor <- update.R2.tau.and.obj(factor, flash)

  if (get.obj(factor) >= get.obj(flash)) {
    flash <- alter.existing.factor(flash, factor)
    flash <- set.nullchk.fail.flag(flash)
  }

  return(flash)
}

zero.factor <- function(flash, k) {
  factor         <- list()
  factor$k       <- k
  factor$EF      <- r1.zeros(flash)
  factor$EF2     <- factor$EF
  factor$KL      <- rep(0, get.dim(flash))
  factor$is.zero <- TRUE

  return(factor)
}
