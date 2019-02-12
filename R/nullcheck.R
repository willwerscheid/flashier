nullcheck.factor <- function(flash, k, verbose.lvl) {
  if (!is.valid(flash, k))
    return(flash)

  factor <- zero.factor(flash, k)
  factor <- update.R2.tau.and.obj(factor, flash)
  obj.diff <- get.obj(factor) - get.obj(flash)

  if (obj.diff >= 0) {
    flash <- alter.existing.factor(flash, factor)
    flash <- set.nullchk.fail.flag(flash)
    report.nullchk.failure(verbose.lvl, obj.diff, k)
  }

  return(flash)
}

zero.factor <- function(flash, k) {
  factor            <- list()
  factor$k          <- k
  factor$EF         <- r1.zeros(flash)
  factor$EF2        <- factor$EF
  factor$g          <- list()
  factor$KL         <- rep(0, get.dim(flash))
  factor$exclusions <- list()
  factor$is.zero    <- TRUE

  return(factor)
}
