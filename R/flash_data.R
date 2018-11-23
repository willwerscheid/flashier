set.flash.data <- function(data, S = NULL, tau = NULL, tau.dim = NULL) {
  if (is(data, "flash.data")) {
    for (arg in list(S, tau, tau.dim)) {
      if (!is.null(arg))
        warning("Ignoring ", deparse(substitute(arg)))
    }
  } else {
    must.not.supply.both(S, tau)
    must.not.supply.both(S, tau.dim)
    must.supply.neither.or.both(tau, tau.dim)
    must.be.supported.data.type(data, allow.null = FALSE)
    must.be.supported.data.type(S, allow.scalar = TRUE)
    must.be.supported.data.type(tau, allow.vector = TRUE)
    must.be.valid.integer(tau.dim, lower = 0, upper = length(dim(data)))

    flash.data        <- list()
    class(flash.data) <- "flash.data"
    flash.data$Y      <- data
    if (anyNA(data)) {
      flash.data$Y[is.na(data)] <- 0
      flash.data$Z              <- 1L * is.na(data)
    } else {
      flash.data$Z <- 1
    }

    if (!is.null(S)) {
      if (length(S) == 1) {
        tau.dim <- 0
      } else {
        dims.must.match(data, S)
      }
      tau <- 1 / S^2
    } else if (!is.null(tau)) {
      dims.must.match(data, tau, tau.dim)
    }
    flash.data$tau     <- tau
    flash.data$tau.dim <- tau.dim
  }

  return(flash.data)
}
