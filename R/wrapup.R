wrapup.flash <- function(flash, output.lvl) {
  class(flash) <- c("flash.fit", "list")

  if (output.lvl == 0) {
    flash <- set.bypass.init.flag(flash)
    return(flash)
  }

  flash.object <- list()

  flash.object$n.factors <- get.n.factors(flash)

  if (flash.object$n.factors > 0) {
    flash.object$pve <- calc.pve(flash)
  }

  flash.object$elbo <- get.obj(flash)

  if (is.tau.simple(flash)) {
    flash.object$residuals.sd <- 1 / sqrt(get.tau(flash))
  } else if (is.list(get.tau(flash))) {
    flash.object$residuals.sd <- lapply(get.tau(flash),
                                        function(tau) 1 / sqrt(tau))
  }

  all.valid <- all(is.valid(flash))

  if (flash.object$n.factors > 0) {
    begin.idx <- length(flash.object) + 1
    if (get.dim(flash) == 2) {
      flash.object$L.pm  <- get.EF(flash)[[1]]
      flash.object$F.pm  <- get.EF(flash)[[2]]
      flash.object$L.psd <- get.EF2(flash)[[1]] - get.EF(flash)[[1]]^2
      flash.object$F.psd <- get.EF2(flash)[[2]] - get.EF(flash)[[2]]^2
      if (all.valid) {
        lfsr <- calc.lfsr(flash)
        if (is.matrix(lfsr[[1]])) {
          flash.object$L.lfsr <- lfsr[[1]]
        }
        if (is.matrix(lfsr[[2]])) {
          flash.object$F.lfsr <- lfsr[[2]]
        }
      }
      for (i in begin.idx:length(flash.object)) {
        if ((i - begin.idx) %% 2 == 0) {
          flash.object[i:(i + 1)] <- propagate.names(flash.object[i:(i + 1)], flash)
        }
      }
    } else {
      flash.object$loadings.pm  <- get.EF(flash)
      flash.object$loadings.psd <- mapply(
        get.EF2(flash), get.EF(flash),
        FUN = function(EX2, EX) {sqrt(pmax(EX2 - EX^2, 0))},
        SIMPLIFY = FALSE
      )
      if (all.valid) {
        flash.object$loadings.lfsr <- calc.lfsr(flash)
      }
      for (i in begin.idx:length(flash.object)) {
        flash.object[[i]] <- propagate.names(flash.object[[i]], flash)
      }
    }
  }

  if (flash.object$n.factors > 0) {
    fitted.g <- get.g.by.mode(flash)

    if (get.dim(flash) == 2) {
      names(fitted.g) <- c("L.ghat", "F.ghat")
      flash.object <- c(flash.object, fitted.g)
    } else {
      flash.object$fitted.g <- fitted.g
    }
  }

  if (flash.object$n.factors > 0 && output.lvl > 1 && all.valid) {
    flash.object$sampler <- build.sampler(flash)
  }

  if (output.lvl < 3) {
    flash <- clear.flags(flash)
    flash <- remove.data.elements(flash)
    flash <- remove.auxiliary.elements(flash)
  }

  flash.object$flash.fit <- flash

  class(flash.object) <- list("flash", "list")

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
      lfsr <- rep(NA, get.dims(flash)[n])
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
