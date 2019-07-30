wrapup.flash <- function(flash, output.lvl, is.converged) {
  class(flash) <- "flash.fit"

  if (output.lvl == 0) {
    flash <- set.bypass.init.flag(flash)
    return(flash)
  }

  flash.object <- list()

  flash.object$n.factors        <- get.n.factors(flash)

  if (flash.object$n.factors > 0) {
    flash.object$pve            <- calc.pve(flash)
    loadings                    <- calc.normalized.loadings(flash)
    flash.object$loadings.scale <- loadings$scale.constants
    flash.object$loadings.pm    <- loadings$normalized.loadings
    flash.object$loadings.psd   <- loadings$loading.SDs
    flash.object$loadings.lfsr  <- calc.lfsr(flash)
  }

  if (is.tau.simple(flash)) {
    flash.object$residuals.sd <- 1 / sqrt(get.tau(flash))
  }

  flash.object$fitted.g <- get.g.by.mode(flash)
  flash.object$elbo     <- get.obj(flash)

  if (is.converged) {
    flash.object$convergence.status <- "converged"
  } else {
    flash.object$convergence.status <- paste("did not converge (reached",
                                             "maximum number of iterations)")
  }

  if (flash.object$n.factors > 0 && output.lvl > 1) {
    flash.object$sampler <- build.sampler(flash)
  }

  if (output.lvl < 3) {
    flash <- clear.flags(flash)
    flash <- remove.data.elements(flash)
    flash <- remove.auxiliary.elements(flash)
  } else {
    flash <- set.bypass.init.flag(flash)
  }

  flash.object$flash.fit <- flash

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
  ldf <- calc.normalized.loadings(flash, for.pve = TRUE)
  S   <- ldf$scale.constants

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

calc.normalized.loadings <- function(flash, for.pve = FALSE) {
  ret <- list()

  if (for.pve) {
    loadings <- get.EF2(flash)
    norms <- lapply(loadings, colSums)
  } else {
    loadings <- get.EF(flash)
    norms <- lapply(loadings, function(x) {sqrt(colSums(x^2))})
  }

  # Zero factors are "normalized" to zero.
  norms <- lapply(norms, function(x) {x[is.zero(flash)] <- Inf; x})
  L <- mapply(loadings, norms, FUN = function(X, y) {
    X / rep(y, each = nrow(X))
  }, SIMPLIFY = FALSE)
  if (!for.pve) {
    L2 <- mapply(get.EF2(flash), norms, FUN = function(X, y) {
      X / rep(y^2, each = nrow(X))
    }, SIMPLIFY = FALSE)
    SD <- mapply(L2, L,
                 FUN = function(EX2, EX) {sqrt(pmax(EX2 - EX^2, 0))},
                 SIMPLIFY = FALSE)
  }

  # Propagate names.
  data.dimnames <- get.dimnames(flash)
  for (n in 1:get.dim(flash)) {
    if (!is.null(data.dimnames) && !is.null(data.dimnames[[n]]))
      rownames(L[[n]]) <- data.dimnames[[n]]
  }

  norms <- do.call(rbind, norms)
  ret$scale.constants <- apply(norms, 2, prod)
  ret$scale.constants[is.zero(flash)] <- 0
  ret$normalized.loadings <- L
  if (!for.pve)
    ret$loading.SDs <- SD

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
    if (!is.null(ebnm.res$posterior) && !is.null(ebnm.res$posterior$lfsr)) {
      lfsr <- ebnm.res$posterior$lfsr
      fix.dim <- get.fix.dim(factor)
      if (!is.null(fix.dim) && (fix.dim == n))
        lfsr[get.fix.idx(factor)] <- NA
    } else {
      lfsr <- NULL
    }
  }
  return(lfsr)
}

get.g.by.mode <- function(flash) {
  g.by.factor <- get.g(flash)
  g.by.mode <- lapply(1:get.dim(flash), FUN = function(n) {
    lapply(g.by.factor, FUN = function(g.k) {
      if (length(g.k) == 0)
        return(NULL)
      else
        return(g.k[[n]])
    })
  })
  return(g.by.mode)
}
