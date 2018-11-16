calc.est.tau <- function(flash, delta.R2 = 0) {
  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

get.tau.r1 <- function(flash, est.tau = NULL) {
  if (is.null(tau))
    est.tau <- get.est.tau(flash)
  n <- get.tau.n(flash)

  tau.r1 <- as.list(rep(1, get.dim(flash)))
  tau.r1[[n]] <- est.tau
  class(tau.r1) <- "r1"

  return(tau.r1)
}

get.tau.lowrank <- function(flash, est.tau = NULL) {
  if (is.null(est.tau))
    est.tau <- get.est.tau(flash)
  n <- get.tau.n(flash)

  tau.lowrank <- lapply(as.list(get.dims(flash)),
                        function(d) {matrix(1, nrow = d, ncol = 1)})
  tau.lowrank[[n]] <- est.tau * tau.lowrank[[n]]
  class(tau.lowrank) <- "lowrank"

  return(tau.lowrank)
}
