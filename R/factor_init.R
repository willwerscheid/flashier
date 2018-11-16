init.factor <- function(flash) {
  factor          <- list()
  factor$EF       <- init.r1(flash)
  factor$EF2      <- r1.square(factor$EF)
  factor$delta.R2 <- calc.delta.R2(factor, flash)
  factor$est.tau  <- calc.est.tau(flash, factor$delta.R2)
  factor$obj      <- -Inf # not yet a valid factor
  factor$is.zero  <- FALSE

  return(factor)
}

init.r1 <- function(flash, tol = 1e-3) {
  r1 <- r1.random(get.dims(flash))
  r1 <- optimize.it(r1,
                    update.fn = update.r1,
                    update.args = list(flash = flash),
                    obj.fn = calc.max.chg.r1,
                    tol = tol)

  return(r1)
}

update.r1 <- function(r1, flash) {
  for (n in 1:length(r1))
    r1 <- update.r1.one.n(r1, n, flash)

  return(r1)
}

update.r1.one.n <- function(r1, n, flash) {
  R  <- get.R(flash)
  Y  <- get.Y(flash)
  Z  <- get.nonmissing(flash)
  EF <- get.EF(flash)

  if (uses.R(flash)) {
    r1[[n]] <- (nmode.prod.r1(R, r1[-n], n)
                / nmode.prod.r1(Z, r1.square(r1[-n]), n))
  } else {
    r1[[n]] <- ((nmode.prod.r1(Y, r1[-n], n)
                 - premult.nmode.prod.r1(Z, EF, r1[-n], n))
                / nmode.prod.r1(Z, r1.square(r1[-n]), n))
  }

  return(r1)
}
