# TODO: handle is.converged better, document, check args
flash.nullcheck <- function(flash,
                            kset = NULL,
                            tol = set.default.tol(flash),
                            verbose.lvl = 1,
                            output.lvl = 3) {
  if (inherits(flash, "flash")) {
    conv.stat <- get.conv.stat(flash)
    flash <- get.fit(flash)
  } else {
    conv.stat <- NULL
  }

  if (is.null(kset)) {
    kset <- 1:get.n.factors(flash)
  }

  announce.nullchk(verbose.lvl, n.factors = length(kset))

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl, tol)
  }

  if (length(kset) > 0 && !nullchk.failed(flash)) {
    report.nullchk.success(verbose.lvl)
  }

  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl, is.converged = TRUE)

  if (nullchk.failed(flash)) {
    flash$convergence.status <- "not converged (nullcheck failed)"
  } else {
    flash$convergence.status <- conv.stat
  }

  report.completion(verbose.lvl)
  return(flash)
}

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
