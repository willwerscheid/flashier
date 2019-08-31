nullcheck.factor <- function(flash, k, verbose.lvl, tol) {
  if (!is.valid(flash, k) || is.zero(flash, k))
    return(flash)

  factor <- zero.factor(flash, k)
  factor <- update.R2.tau.and.obj(factor, flash)
  obj.diff <- get.obj(factor) - get.obj(flash)

  if (obj.diff >= -tol) {
    flash <- insert.factor(flash, factor)
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

#' @export
remove.factors <- function(fl, kset, output.lvl = 3) {
  if (inherits(fl, "flash")) {
    fit <- get.fit(fl)
  } else if (inherits(fl, "flash.fit")) {
    fit <- fl
  } else {
    stop("Parameter fl must be a flash or flash.fit object.")
  }

  for (k in kset) {
    fit <- nullcheck.factor(fit, k, verbose.lvl = 0, tol = Inf)
  }

  flash <- wrapup.flash(flash, output.lvl, is.converged = TRUE)

  return(flash)
}
