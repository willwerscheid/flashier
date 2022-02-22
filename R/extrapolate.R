default.extrapolate.param <- list(beta.init = 0.5,
                                  beta.increase = 1.2,
                                  beta.decrease = 0.5,
                                  beta.max = 2)

#' @importFrom utils modifyList
set.extrapolate.param <- function(control) {
  par <- default.extrapolate.param
  if (!all(names(control) %in% names(par))) {
    stop("Unrecognized extrapolation control parameter.")
  }
  par <- modifyList(par, control)
  return(par)
}

extrapolate <- function(new, old, beta) {
  return(new + beta * (new - old))
}

# Works for both factor and flash objects.
extrapolate.f <- function(f, old.f, par) {
  beta <- par$beta
  epsilon <- 1e-10

  EF  <- mapply(extrapolate, get.EF(f), get.EF(old.f),
                MoreArgs = list(beta = beta), SIMPLIFY = FALSE)
  EF2 <- mapply(extrapolate, get.EF2(f), get.EF2(old.f),
                MoreArgs = list(beta = beta), SIMPLIFY = FALSE)

  # Ensure that EF2 > EF^2.
  EF2 <- mapply(function(EF, EF2) pmax(EF2, EF^2 + epsilon), EF, EF2,
                SIMPLIFY = FALSE)

  if (is.list(get.tau(f))) {
    tau <- mapply(extrapolate, get.tau(f), get.tau(old.f),
                  MoreArgs = list(beta = beta), SIMPLIFY = FALSE)
    # Ensure that tau > 0.
    tau <- lapply(tau, function(tau) pmax(tau, epsilon))
  } else {
    tau <- extrapolate(get.tau(f), get.tau(old.f), beta)
    tau <- pmax(tau, epsilon)
  }

  if (!is.tau.lowrank(f)) {
    tau <- tau * get.nonmissing(f)
  }

  class(EF)  <- class(get.EF(f))
  class(EF2) <- class(get.EF2(f))
  class(tau) <- class(get.tau(f))

  f <- set.EF(f, EF)
  f <- set.EF2(f, EF2)
  f <- set.tau(f, tau)

  if (uses.R(f))
    f <- set.R(f, get.Y(f) - lowrank.expand(get.EF(f)))

  return(f)
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
