#' @export
#
flash.remove.factors <- function(flash, kset, output.lvl = 3) {
  if (inherits(flash, "flash")) {
    conv.stat <- get.conv.stat(flash)
    flash <- get.fit(flash)
  } else {
    conv.stat <- NULL
  }

  for (k in kset) {
    flash <- nullcheck.factor(flash, k, verbose.lvl = 0, tol = Inf)
  }

  flash <- wrapup.flash(flash, output.lvl, is.converged = TRUE)

  return(flash)
}
