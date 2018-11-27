wrapup.flash <- function(flash, data) {
  flash.object <- list()
  flash.object$n.factors <- get.n.factors(flash)
  flash.object$pve       <- calc.pve(flash)
  flash.object$loadings  <- calc.normalized.loadings(flash)
  flash.object$objective <- get.obj(flash)
  flash.object$sampler   <- F.sampler(flash)

  flash <- clear.flags(flash)
  flash <- remove.data.elements(flash)
  flash <- remove.auxiliary.elements(flash)
  class(flash) <- "flash.fit"

  flash.object$fit <- flash
  class(flash.object) <- "flash"
  return(flash.object)
}

remove.data.elements <- function(flash) {
  flash <- set.R(flash, NULL)
  flash <- set.Y(flash, NULL)
  flash <- set.nonmissing(flash, NULL)
  return(flash)
}

# TODO: remove est.tau once tests no longer need it
remove.auxiliary.elements <- function(flash) {
  flash <- set.n.nonmissing(flash, NULL)
  flash <- set.R2(flash, NULL)
  # flash <- set.est.tau(flash, NULL)
  return(flash)
}

calc.pve <- function(flash) {
  ldf <- calc.normalized.loadings(flash)
  S   <- ldf$scale.factor^2

  tau <- get.tau(flash)
  if (is.tau.simple(flash)) {
    var.from.tau <- sum(get.n.nonmissing(flash) / tau)
  } else if (is.var.type.kronecker(flash)) {
    var.from.tau <- sum(get.nonmissing(flash) / r1.expand(tau))
  } else{
    var.from.tau <- sum(1 / tau[tau > 0])
  }

  return(S / (sum(S) + var.from.tau))
}

calc.normalized.loadings <- function(flash) {
  ret <- list()

  EF <- get.EF(flash)
  norms <- lapply(EF, function(x) {sqrt(colSums(x^2))})
  L <- mapply(EF, norms, FUN = function(X, y) {
    X / matrix(y, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  })

  if (uses.R(flash)) {
    data.dimnames <- dimnames(get.R(flash))
  } else {
    data.dimnames <- dimnames(get.Y(flash))
  }
  if (!is.null(data.dimnames)) {
    for (n in 1:get.dim(flash)) {
      rownames(L[[n]]) <- data.dimnames[[n]]
    }
  }

  norms <- do.call(rbind, norms)
  ret$scale.factor <- apply(norms, 2, prod)
  ret$normalized.loadings <- L

  return(ret)
}

