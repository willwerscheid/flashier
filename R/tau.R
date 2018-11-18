calc.est.tau <- function(flash, delta.R2 = 0) {
  return(get.n.nonmissing(flash) / (get.R2(flash) + delta.R2))
}

get.tau.lowrank <- function(flash, est.tau = NULL) {
  if (is.null(est.tau))
    est.tau <- get.est.tau(flash)
  n <- get.tau.n(flash)

  tau.lowrank <- lapply(as.list(get.dims(flash)),
                        function(dim) {matrix(1, nrow = dim, ncol = 1)})
  tau.lowrank[[n]] <- est.tau * tau.lowrank[[n]]
  class(tau.lowrank) <- "lowrank"

  return(tau.lowrank)
}
