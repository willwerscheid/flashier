init.factor <- function(flash, init.fn) {
  factor <- list()

  factor$EF <- do.call(init.fn, list(flash))
  class(factor$EF) <- c("r1", "list")

  factor$EF2        <- r1.square(factor$EF)
  factor$KL         <- rep(0, get.dim(flash))
  factor            <- update.tau(factor, flash)
  factor$obj        <- calc.obj(flash, factor)
  factor$is.zero    <- all(unlist(factor$EF) == 0)
  factor$is.valid   <- factor$is.zero

  return(factor)
}
