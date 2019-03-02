# TODO: Change default to output LFSR after ebnm can do it and after
#   changing default ash method to "fdr".

wrapup.flash <- function(flash, output.lvl) {
  class(flash) <- "flash.fit"

  if (output.lvl == 0) {
    flash <- set.bypass.init.flag(flash)
    return(flash)
  }

  flash.object <- list()

  flash.object$n.factors   <- get.n.factors(flash)
  flash.object$objective   <- get.obj(flash)
  if (flash.object$n.factors > 0) {
    flash.object$pve       <- calc.pve(flash)
    flash.object$loadings  <- calc.normalized.loadings(flash)
    if (output.lvl > 3)
      flash.object$lfsr    <- calc.lfsr(flash)
    if (output.lvl > 1)
      flash.object$sampler <- F.sampler(flash)
  }

  if (output.lvl < 3) {
    flash <- clear.flags(flash)
    flash <- remove.data.elements(flash)
    flash <- remove.auxiliary.elements(flash)
  } else {
    flash <- set.bypass.init.flag(flash)
  }

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

remove.auxiliary.elements <- function(flash) {
  flash <- set.n.nonmissing(flash, NULL)
  flash <- set.R2(flash, NULL)
  flash <- set.est.tau(flash, NULL)

  return(flash)
}

calc.pve <- function(flash) {
  ldf <- calc.normalized.loadings(flash)
  S   <- ldf$scale.constant^2

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
  # Zero factors are "normalized" to zero.
  norms <- lapply(norms, function(x) {x[is.zero(flash)] <- Inf; x})
  L <- mapply(EF, norms, FUN = function(X, y) {
    X / rep(y, each = nrow(X))
  })

  # Propagate names.
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
  ret$scale.constant <- apply(norms, 2, prod)
  ret$scale.constant[is.zero(flash)] <- 0
  ret$normalized.loadings <- L

  return(ret)
}

calc.lfsr <- function(flash) {
  return(lapply(1:get.dim(flash),
                function(n) sapply(1:get.n.factors(flash),
                                   function(k) lfsr.one.n(flash, k, n))))

}

lfsr.one.n <- function(flash, k, n) {
  factor <- extract.factor(flash, k)
  if (is.zero(factor) || all.fixed(factor, n)) {
    lfsr <- rep(NA, get.dims(flash)[n])
  } else {
    ebnm.res <- solve.ebnm(factor, n, flash, output = "lfsr")
    lfsr <- ebnm.res$lfsr
    if (get.fix.dim(factor) == n)
      lfsr[get.fix.idx(factor)] <- NA
  }
  return(lfsr)
}
