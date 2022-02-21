# Helpers for the flash (not flash.fit) object --------------------------------

get.fit <- function(f) {
  if (inherits(f, "flash"))
    return(f[["flash.fit"]])
  if (inherits(f, "flash.fit"))
    return(f)
  stop("flash must be a flash or flash.fit object. Use flash.init to ",
       "initialize a flash.fit object.")
}
set.fit <- function(f, fit) {
  f[["flash.fit"]] <- fit
  return(f)
}
get.conv.stat <- function(f) f[["convergence.status"]]


# Getters for the flash.fit object (also used by the smaller factors) ---------

get.R                 <- function(f) f[["R"]]
get.nonmissing        <- function(f) f[["Z"]]
get.given.S2          <- function(f) f[["given.S2"]]
get.given.tau         <- function(f) f[["given.tau"]]
get.given.tau.dim     <- function(f) f[["given.tau.dim"]]
get.est.S2            <- function(f) f[["est.S2"]]
get.est.tau           <- function(f) f[["est.tau"]]
get.est.tau.dim       <- function(f) f[["est.tau.dim"]]
get.tau               <- function(f) f[["tau"]]
use.fixed.to.est.g    <- function(f) f[["use.fixed.to.est.g"]]
get.n.nonmissing      <- function(f) f[["n.nonmissing"]]
get.kron.nonmissing   <- function(f) f[["kron.nonmissing"]]
get.Y2                <- function(f) f[["Y2"]]
get.R2                <- function(f) f[["R2"]]
get.delta.R2          <- function(f) f[["delta.R2"]]
get.log.2pi.s2        <- function(f) f[["log.2pi.s2"]]
get.sum.tau.R2        <- function(f) f[["sum.tau.R2"]]
get.obj               <- function(f) f[["obj"]]
warmstart.backfits    <- function(f) f[["warmstart.backfits"]]
warmstart.greedy      <- function(f) f[["warmstart.greedy"]]

get.Y <- function(f, require.fullrank = FALSE) {
  Y <- f[["Y"]]
  if (require.fullrank && inherits(Y, "lowrank"))
    Y <- lowrank.expand(Y)
  return(Y)
}
get.EF <- function(f, n = NULL) {
  EF <- f[["EF"]]
  if (is.null(EF[[1]]))
    return(NULL)
  if (!is.null(n))
    return(EF[[n]])
  return(EF)
}
get.EF.k <- function(f, k, n = NULL) {
  EFk <- lapply(f[["EF"]], function(X) X[, k])
  if (!is.null(n))
    return(EFk[[n]])
  class(EFk) <- "r1"
  return(EFk)
}
get.EF2 <- function(f, n = NULL) {
  EF2 <- f[["EF2"]]
  if (is.null(EF2[[1]]))
    return(NULL)
  if (!is.null(n))
    return(EF2[[n]])
  return(EF2)
}
get.EF2.k <- function(f, k) {
  EF2k <- lapply(f[["EF2"]], function(X) X[, k])
  class(EF2k) <- "r1"
  return(EF2k)
}
get.dim.signs <- function(f, k = NULL) {
  if (is.null(k))
    return(unlist(f[["dim.signs"]]))
  return(unlist(f[["dim.signs"]][[min(k, length(f[["dim.signs"]]))]]))
}
get.fix.dim <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["fix.dim"]])
  if (length(f[["fix.dim"]]) < k)
    return(NULL)
  return(f[["fix.dim"]][[k]])
}
get.fix.idx <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["fix.idx"]])
  if (length(f[["fix.idx"]]) < k)
    return(NULL)
  return(f[["fix.idx"]][[k]])
}
get.fix.vals <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["fix.vals"]])
  if (length(f[["fix.vals"]]) < k)
    return(NULL)
  return(f[["fix.vals"]][[k]])
}
get.ebnm.fn.k <- function(f, k) {
  if (!is.list(f[["ebnm.fn"]]))
    return(f[["ebnm.fn"]])
  return(f[["ebnm.fn"]][[k]])
}
get.g <- function(f, n = NULL) {
  if (is.null(n))
    return(f[["g"]])
  if (is.null(f[["g"]]) || length(f[["g"]]) < n)
    return(NULL)
  return(f[["g"]][[n]])
}
get.g.k <- function(f, k, n = NULL) {
  if (is.null(k) && is.null(n))
    return(f[["g"]])
  if (is.null(k))
    return(lapply(f[["g"]], `[[`, n))
  if (is.null(n))
    return(f[["g"]][[k]])
  return(f[["g"]][[k]][[n]])
}
get.KL <- function(f, n = NULL) {
  if (is.null(n))
    return(f[["KL"]])
  if (is.null(f[["KL"]]))
    return(NULL)
  return(f[["KL"]][[n]])
}
get.KL.k <- function(f, k, n = NULL) {
  if (is.null(k)) {
    KL <- f[["KL"]]
  } else {
    KL <- sapply(f[["KL"]], getElement, k)
  }
  if (is.null(n))
    return(KL)
  return(KL[[n]])
}
is.zero <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["is.zero"]])
  return(f[["is.zero"]][k])
}
is.valid <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["is.valid"]])
  return(f[["is.valid"]][k])
}
get.nonmissing.thresh <- function(f, n) f[["nonmissing.thresh"]][n]
get.exclusions <- function(f, n = NULL) {
  if (is.null(n))
    return(f[["exclusions"]])
  if (length(f[["exclusions"]]) < n)
    return(NULL)
  return(f[["exclusions"]][[n]])
}
get.verbose.lvl <- function(f) {
  f <- get.fit(f)
  f[["verbose.lvl"]]
}
get.verbose.fns <- function(f) {
  f <- get.fit(f)
  f[["verbose.fns"]]
}
get.verbose.colnames <- function(f) {
  f <- get.fit(f)
  f[["verbose.colnames"]]
}
get.verbose.colwidths <- function(f) {
  f <- get.fit(f)
  f[["verbose.colwidths"]]
}


# Additional getters that are only used by factors ----------------------------

get.k          <- function(f) f[["k"]]
get.R.subset   <- function(f) f[["subset.data"]][["R.subset"]]
get.Y.subset   <- function(f) f[["subset.data"]][["Y.subset"]]
get.Z.subset   <- function(f) f[["subset.data"]][["Z.subset"]]
get.EF.subset  <- function(f) f[["subset.data"]][["EF.subset"]]
get.EF2.subset <- function(f) f[["subset.data"]][["EF2.subset"]]

get.idx.subset   <- function(f) {
  if (!is.null(f[["subset.data"]]))
    return(f[["subset.data"]][["idx.subset"]])
  return(f[["idx.subset"]])
}
is.new <- function(f) is.null(get.k(f))


# Simple helper functions for the main flash object and smaller factors -------

get.n.factors <- function(f) {
  max(0, ncol(f[["EF"]][[1]]))
}
get.dims <- function(f) {
  if (uses.R(f))
    return(dim(get.R(f)))
  Y <- get.Y(f)
  if (inherits(Y, "lowrank"))
    return(sapply(Y, nrow))
  return(dim(Y))
}
get.dim <- function(f) length(get.dims(f))
get.dimnames <- function(f) {
  if (uses.R(f))
    return(dimnames(get.R(f)))
  Y <- get.Y(f)
  if (inherits(Y, "lowrank"))
    return(lapply(Y, rownames))
  return(dimnames(Y))
}
get.next.k <- function(f) {
  return(get.n.factors(f) + 1)
}
any.missing <- function(f) !identical(get.nonmissing(f), 1)
is.obj.valid <- function(flash, factor = NULL) {
  valid <- is.valid(flash)
  if (!is.null(factor))
    valid <- c(valid, is.valid(factor))
  return(all(valid))
}

get.new.EF <- function(flash, factor = NULL) {
  EF <- get.EF(flash)
  if (is.null(factor))
    return(EF)
  if (!is.new(factor))
    EF <- lowrank.drop.k(EF, get.k(factor))
  return(lowranks.combine(EF, as.lowrank(get.EF(factor))))
}
get.new.EF2 <- function(flash, factor = NULL) {
  EF2 <- get.EF2(flash)
  if (is.null(factor))
    return(EF2)
  if (!is.new(factor))
    EF2 <- lowrank.drop.k(EF2, get.k(factor))
  return(lowranks.combine(EF2, as.lowrank(get.EF2(factor))))
}

get.ebnm.fn <- function(flash, factor = NULL, n = NULL) {
  if (is.null(factor))
    return(flash[["ebnm.fn"]])
  if (is.new(factor)) {
    ebnm.fn <- factor[["ebnm.fn"]]
  } else {
    ebnm.fn <- get.ebnm.fn.k(flash, get.k(factor))
  }
  if (is.null(n))
    return(ebnm.fn)
  if (length(ebnm.fn) == 1)
    return(ebnm.fn)
  return(ebnm.fn[[n]])
}
extend.ebnm.lists <- function(flash) {
  k <- get.n.factors(flash)
  l <- length(flash[["ebnm.fn"]])
  extend.by <- (k + 1) - l
  if (extend.by > 0) {
    flash[["ebnm.fn"]] <- c(flash[["ebnm.fn"]],
                            rep(flash[["ebnm.fn"]][l], extend.by))
  }
  return(flash)
}

is.var.type.zero <- function(f) {
  return(is.null(get.est.tau.dim(f)))
}
is.var.type.kronecker <- function(f) {
  return(is.null(get.given.tau(f))
         && is.null(get.given.S2(f))
         && (length(get.est.tau.dim(f)) > 1))
}
is.var.type.noisy <- function(f) {
  return(!is.null(get.given.S2(f))
         && (length(get.est.tau.dim(f)) == 1))
}
is.var.type.noisy.kron <- function(f) {
  return(!is.null(get.given.S2(f))
         && (length(get.est.tau.dim(f)) > 1))
}
is.tau.constant <- function(f) {
  return(!is.null(get.est.tau.dim(f)) && (get.est.tau.dim(f) == 0))
}
is.tau.simple <- function(f) {
  if (is.var.type.noisy(f) || is.var.type.noisy.kron(f))
    return(FALSE)
  var.type <- get.est.tau.dim(f)
  SEs <- get.given.tau(f)
  SEs.dim <- get.given.tau.dim(f)
  # Zero var.type with tau representable as low-rank:
  is.simple <- (is.var.type.zero(f) && !is.null(SEs.dim))
  # Simple by.row/by.column estimation with SEs not provided:
  is.simple <- (is.simple || ((length(var.type) == 1) && is.null(SEs)))
  # by.row/by.column var.type where SEs lie in same dimension as estimated var:
  is.simple <- (is.simple
                || ((length(var.type) == 1)
                    && !is.null(SEs.dim)
                    && (SEs.dim %in% c(0, var.type))))
  return(is.simple)
}
is.tau.lowrank <- function(f) {
  tau <- get.given.tau(f)
  if (is.null(tau))
    tau <- get.given.S2(f)
  return(is.null(tau) || is.vector(tau))
}

get.R2.n <- function(f) {
  n <- max(get.est.tau.dim(f), 0)
  if (n == 0)
    n <- which.min(get.dims(f))
  return(n)
}
store.R2.as.scalar <- function(f) is.tau.constant(f)
store.R2.as.matrix <- function(f) is.var.type.kronecker(f)
uses.R <- function(f) !is.null(get.R(f))
get.new.Rsquared <- function(flash, factor = NULL, EF = NULL,
                             set.missing.to.zero = TRUE) {
  if (uses.R(flash) && !is.null(factor)) {
    R2 <- get.R(factor)^2
  } else if (uses.R(flash)) {
    R2 <- get.R(flash)^2
  } else {
    if (is.null(EF))
      EF <- get.new.EF(flash, factor)
    R2 <- (get.Y(flash, require.fullrank = TRUE) - lowrank.expand(EF))^2
    if (set.missing.to.zero)
      R2 <- get.nonmissing(flash) * R2
  }
  return(R2)
}

get.n.fixed.to.add <- function(f) {
  return(sum(which.k.fixed(f) > get.n.factors(f)))
}
which.k.fixed <- function(f) {
  fix.dim <- get.fix.dim(f)
  if (is.null(fix.dim))
    return(NULL)
  not.fixed <- sapply(fix.dim, is.null)
  return(which(!not.fixed))
}
is.next.fixed <- function(f) {
  return(!is.null(get.fix.dim(f, get.next.k(f))))
}
get.unfixed.idx <- function(f, k) {
  fix.dim <- get.fix.dim(f, k)
  fix.idx <- get.fix.idx(f, k)
  return(setdiff(1:(get.dims(f)[[fix.dim]]), fix.idx))
}
get.next.unfixed.idx <- function(f) {
  return(get.unfixed.idx(f, get.next.k(f)))
}
all.fixed <- function(f, n) {
  fix.dim <- as.integer(get.fix.dim(f))
  idx.subset <- get.idx.subset(f)
  return(identical(fix.dim, n) && length(idx.subset) == 0)
}

greedy.failed <- function(f) {
  return(identical(f[["greedy.fail"]], TRUE))
}
nullchk.failed <- function(f) {
  return(identical(f[["nullchk.fail"]], TRUE))
}
bypass.init <- function(f) {
  return(identical(f[["bypass.init"]], TRUE))
}

get.subset.data <- function(f, fix.dim, idx.subset) {
  if (length(idx.subset) < 1)
    return(NULL)
  subset.data <- list(idx.subset = idx.subset)
  subset.data$R.subset  <- fullrank.subset(get.R(f), fix.dim, idx.subset)
  subset.data$Y.subset  <- full.or.lowrank.subset(get.Y(f), fix.dim, idx.subset)
  subset.data$Z.subset  <- fullrank.subset(get.nonmissing(f), fix.dim, idx.subset)
  subset.data$EF.subset <- lowrank.subset(get.EF(f), fix.dim, idx.subset)
  return(subset.data)
}


# Setters for the main flash object and smaller factors -----------------------

set.R <- function(f, R) {
  f[["R"]] <- R
  return(f)
}
set.Y <- function(f, Y) {
  f[["Y"]] <- Y
  return(f)
}
set.nonmissing <- function(f, Z) {
  f[["Z"]] <- Z
  return(f)
}
set.EF <- function(f, EF, n = NULL) {
  if (is.null(n)) {
    f[["EF"]] <- EF
  } else {
    f[["EF"]][[n]] <- EF
  }
  return(f)
}
set.EFk <- function(f, k, EF) {
  for (n in 1:length(EF)) {
    f[["EF"]][[n]][, k] <- EF[[n]]
  }
  return(f)
}
add.factor.to.EF <- function(f, new.EF) {
  if (is.null(f[["EF"]])) {
    f[["EF"]] <- as.lowrank(new.EF)
  } else {
    f[["EF"]] <- lowranks.combine(f[["EF"]], new.EF)
  }
  return(f)
}
set.EF2 <- function(f, EF2, n = NULL) {
  if (is.null(n)) {
    f[["EF2"]] <- EF2
  } else {
    f[["EF2"]][[n]] <- EF2
  }
  return(f)
}
set.EF2k <- function(f, k, EF2) {
  for (n in 1:length(EF2)) {
    f[["EF2"]][[n]][, k] <- EF2[[n]]
  }
  return(f)
}
add.factor.to.EF2 <- function(f, new.EF2) {
  if (is.null(f[["EF2"]])) {
    f[["EF2"]] <- as.lowrank(new.EF2)
  } else {
    f[["EF2"]] <- lowranks.combine(f[["EF2"]], new.EF2)
  }
  return(f)
}
set.KL <- function(f, KL, n = NULL) {
  if (is.null(n)) {
    f[["KL"]] <- KL
  } else {
    if (is.null(f[["KL"]]))
      f[["KL"]] <- rep(list(NULL), get.dim(f))
    f[["KL"]][[n]] <- KL
  }
  return(f)
}
set.KLk <- function(f, k, KL) {
  for (n in 1:length(KL)) {
    f[["KL"]][[n]][k] <- KL[[n]]
  }
  return(f)
}
add.factor.to.KL <- function(f, new.KL) {
  if (is.null(f[["KL"]])) {
    f[["KL"]] <- as.list(new.KL)
  } else {
    f[["KL"]] <- mapply(c, get.KL(f), new.KL, SIMPLIFY = FALSE)
  }
  return(f)
}
set.g <- function(f, g, n = NULL) {
  if (is.null(n)) {
    f[["g"]] <- g
  } else {
    f[["g"]][[n]] <- g
  }
  return(f)
}
set.gk <- function(f, k, g) {
  f[["g"]][[k]] <- g
  return(f)
}
add.factor.to.g <- function(f, new.g) {
  f[["g"]] <- c(f[["g"]], list(new.g))
  return(f)
}
set.ebnm.fn <- function(f, ebnm.fn) {
  f[["ebnm.fn"]] <- ebnm.fn
  return(f)
}
add.factor.to.ebnm.fn <- function(f, new.ebnm.fn) {
  f[["ebnm.fn"]] <- c(f[["ebnm.fn"]], list(new.ebnm.fn))
  return(f)
}
set.est.S2 <- function(f, est.S2) {
  f[["est.S2"]] <- est.S2
  return(f)
}
set.est.tau <- function(f, est.tau) {
  f[["est.tau"]] <- est.tau
  return(f)
}
set.tau <- function(f, tau) {
  f[["tau"]] <- tau
  return(f)
}
set.obj <- function(f, obj) {
  f[["obj"]] <- obj
  return(f)
}
set.n.nonmissing <- function(f, n.nonmissing) {
  f[["n.nonmissing"]] <- n.nonmissing
  return(f)
}
set.Y2 <- function(f, Y2) {
  f[["Y2"]] <- Y2
  return(f)
}
set.R2 <- function(f, R2) {
  f[["R2"]] <- R2
  return(f)
}
set.delta.R2 <- function(f, delta.R2) {
  f[["delta.R2"]] <- delta.R2
  return(f)
}
set.sum.tau.R2 <- function(f, sum.tau.R2) {
  f[["sum.tau.R2"]] <- sum.tau.R2
  return(f)
}
set.fixed.to.est.g <- function(f, use.fixed) {
  f[["use.fixed.to.est.g"]] <- use.fixed
}
set.is.zero <- function(f, is.zero) {
  f[["is.zero"]] <- is.zero
  return(f)
}
add.is.zero <- function(f, is.zero) {
  f[["is.zero"]] <- c(f[["is.zero"]], is.zero)
  return(f)
}
set.to.zero <- function(f, k = NULL) {
  if (is.null(k)) {
    f[["is.zero"]] <- TRUE
  } else {
    f[["is.zero"]][k] <- TRUE
  }
  return(f)
}
set.is.valid <- function(f, is.valid) {
  f[["is.valid"]] <- is.valid
  return(f)
}
add.is.valid <- function(f, is.valid) {
  f[["is.valid"]] <- c(f[["is.valid"]], is.valid)
  return(f)
}
set.to.valid <- function(f, k = NULL) {
  if (is.null(k)) {
    f[["is.valid"]] <- TRUE
  } else {
    f[["is.valid"]][k] <- TRUE
  }
  return(f)
}
set.fix.dim <- function(f, fix.dim) {
  f[["fix.dim"]] <- fix.dim
  return(f)
}
set.fix.idx <- function(f, fix.idx) {
  f[["fix.idx"]] <- fix.idx
  return(f)
}
set.exclusions <- function(f, exclusions, n) {
  if (is.null(get.k(f))) {
    f[["exclusions"]][[n]] <- exclusions
  } else {
    # During backfitting, exclusions can be added but not removed.
    f[["exclusions"]][[n]] <- union(f[["exclusions"]][[n]], exclusions)
  }
  return(f)
}
add.exclusions <- function(f, exclusions, k = NULL) {
  if (is.null(k)) {
    f[["exclusions"]] <- c(f[["exclusions"]], list(exclusions))
  } else {
    f[["exclusions"]][[k]] <- exclusions
  }
  return(f)
}
set.greedy.fail.flag <- function(f) {
  f[["greedy.fail"]] <- TRUE
  return(f)
}
clear.greedy.fail.flag <- function(f) {
  f[["greedy.fail"]] <- NULL
  return(f)
}
set.nullchk.fail.flag <- function(f) {
  f[["nullchk.fail"]] <- TRUE
  return(f)
}
clear.nullchk.fail.flag <- function(f) {
  f[["nullchk.fail"]] <- NULL
  return(f)
}
clear.flags <- function(f) {
  f <- clear.greedy.fail.flag(f)
  f <- clear.nullchk.fail.flag(f)
  return(f)
}
set.bypass.init.flag <- function(f) {
  f[["bypass.init"]] <- TRUE
  return(f)
}
clear.bypass.init.flag <- function(f) {
  f[["bypass.init"]] <- NULL
  return(f)
}
add.subset.data <- function(factor, flash, fix.dim, idx.subset) {
  factor[["subset.data"]] <- get.subset.data(flash, fix.dim, idx.subset)
  factor[["idx.subset"]]  <- NULL
  return(factor)
}
set.warmstart <- function(f, warmstart) {
  f[["warmstart.backfits"]] <- warmstart
  return(f)
}
set.gwarmstart <- function(f, warmstart) {
  f[["warmstart.greedy"]] <- warmstart
  return(f)
}
set.verbose.options <- function(f, lvl, fns, colnames, colwidths) {
  f[["verbose.lvl"]] <- lvl
  f[["verbose.fns"]] <- fns
  f[["verbose.colnames"]] <- colnames
  f[["verbose.colwidths"]] <- colwidths
  return(f)
}
set.max.backfit.iter.reached.flag <- function(f) {
  f[["maxiter.reached"]] <- TRUE
  return(f)
}
clear.max.backfit.iter.reached.flag <- function(f) {
  f[["maxiter.reached"]] <- NULL
  return(f)
}

# Testing function that converts a flashier object into a flashr fit object.
to.flashr <- function(f) {
  if (inherits(f, "flash"))
    f <- get.fit(f)

  flash      <- list()
  flash$EL   <- f$EF[[1]]
  flash$EF   <- f$EF[[2]]
  flash$EL2  <- f$EF2[[1]]
  flash$EF2  <- f$EF2[[2]]
  flash$fixl <- matrix(FALSE, nrow = nrow(flash$EL), ncol = ncol(flash$EL))
  flash$fixf <- matrix(FALSE, nrow = nrow(flash$EF), ncol = ncol(flash$EF))
  flash$gl   <- lapply(f$g, function(k) k[[1]])
  flash$gf   <- lapply(f$g, function(k) k[[2]])
  flash$KL_l <- as.list(f$KL[[1]])
  flash$KL_f <- as.list(f$KL[[2]])
  flash$tau  <- f$tau

  if (!is.null(f$fix.dim)) {
    for (k in 1:get.n.factors(f)) {
      if (f$fix.dim[[k]] == 1) {
        flash$fixl[f$fix.idx[[k]], k] <- TRUE
      }
      if (f$fix.dim[[k]] == 2) {
        flash$fixf[f$fix.idx[[k]], k] <- TRUE
      }
    }
  }

  class(flash) <- "flash_fit"

  return(flash)
}
