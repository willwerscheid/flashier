extrapolate <- function(new, old, beta) {
  return(new + beta * (new - old))
}

extrapolate.factor <- function(factor, old.f, par) {
  beta <- par$beta
  epsilon <- 1e-10

  EF  <- mapply(extrapolate, get.EF(factor), get.EF(old.f),
                MoreArgs = list(beta = beta), SIMPLIFY = FALSE)
  EF2 <- mapply(extrapolate, get.EF2(factor), get.EF2(old.f),
                MoreArgs = list(beta = beta), SIMPLIFY = FALSE)

  # Ensure that EF2 > EF^2.
  EF2 <- mapply(function(EF, EF2) pmax(EF2, EF^2 + epsilon), EF, EF2,
                SIMPLIFY = FALSE)

  class(EF)  <- "r1"
  class(EF2) <- "r1"

  if (is.list(get.tau(factor))) {
    tau <- mapply(extrapolate, get.tau(factor), get.tau(old.f),
                  MoreArgs = list(beta = beta), SIMPLIFY = FALSE)
    tau <- lapply(tau, function(tau) pmax(tau, epsilon))
    class(tau) <- "r1"
  } else {
    tau <- extrapolate(get.tau(factor), get.tau(old.f), beta)
    tau <- pmax(tau, epsilon)
  }

  factor <- set.EF(factor, EF)
  factor <- set.EF2(factor, EF2)
  factor <- set.tau(factor, tau)

  return(factor)
}

init.beta <- function(par) {
  # First iteration always "succeeds", so start slower.
  par$beta <- par$beta.init / par$beta.increase
  return(par)
}

accelerate <- function(par) {
  par$beta <- min(par$beta * par$beta.increase, par$beta.max)
  return(par)
}

decelerate <- function(par) {
  par$beta <- par$beta * par$beta.decrease
  return(par)
}

default.extrapolate.param <- function() {
  return(list(beta.init = 0.5,
              beta.increase = 1.2,
              beta.decrease = 0.5,
              beta.max = 2))
}

set.extrapolate.param <- function(control) {
  par <- default.extrapolate.param()
  if (!all(names(control) %in% names(par))) {
    stop("Unrecognized extrapolation control parameter.")
  }
  par <- modifyList(par, control)
  return(par)
}
