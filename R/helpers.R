get.R            <- function(f) f[["R"]]
get.Y            <- function(f) f[["Y"]]
get.nonmissing   <- function(f) f[["Z"]]
get.g            <- function(f) f[["g"]]
get.KL           <- function(f) f[["KL"]]
get.tau.dim      <- function(f) f[["tau.dim"]]
get.tau.n        <- function(f) f[["tau.n"]]
get.est.tau      <- function(f) f[["est.tau"]]
get.R2           <- function(f) f[["R2"]]
get.n.nonmissing <- function(f) f[["n.nonmissing"]]
get.obj          <- function(f) f[["obj"]]
get.ebnm.fn      <- function(f) f[["ebnm.fn"]]
get.ebnm.param   <- function(f) f[["ebnm.param"]]
get.k            <- function(f) f[["k"]]
get.delta.R2     <- function(f) f[["delta.R2"]]
get.R.subset     <- function(f) f[["subset.data"]][["R.subset"]]
get.Y.subset     <- function(f) f[["subset.data"]][["Y.subset"]]
get.Z.subset     <- function(f) f[["subset.data"]][["Z.subset"]]
get.EF.subset    <- function(f) f[["subset.data"]][["EF.subset"]]
get.EF2.subset   <- function(f) f[["subset.data"]][["EF2.subset"]]
is.fixed         <- function(f) f[["is.fixed"]]
is.valid         <- function(f) f[["is.valid"]]

is.new             <- function(f) is.null(f[["k"]])
store.R2.as.scalar <- function(f) get.tau.dim(f) < 1
use.fixed.to.est.g <- function(f) f[["use.fixed.to.est.g"]]

get.EF <- function(f, n = NULL) {
  EF <- f[["EF"]]
  if (is.null(EF[[1]]))
    return(NULL)
  if (!is.null(n))
    return(EF[[n]])
  return(EF)
}
get.EF2 <- function(f, n = NULL) {
  EF2 <- f[["EF2"]]
  if (is.null(EF2[[1]]))
    return(NULL)
  if (!is.null(n))
    return(EF2[[n]])
  return(EF2)
}

uses.R <- function(f) !is.null(f[["R"]])

is.zero <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["is.zero"]])
  return(f[["is.zero"]][k])
}

is.obj.valid <- function(flash, factor = NULL) {
  valid <- c(flash[["is.valid"]])
  if (!is.null(factor))
    valid <- c(valid, factor[["is.valid"]])
  return(all(valid))
}

get.n.factors <- function(f) max(0, ncol(f$EF[[1]]))
get.dims <- function(f) {
  if (!is.null(get.R(f)))
    return(dim(get.R(f)))
  return(dim(get.Y(f)))
}
get.dim <- function(f) length(get.dims(f))

get.EFk <- function(f, k) {
  EFk <- lapply(f[["EF"]], function(X) X[, k])
  class(EFk) <- "r1"
  return(EFk)
}
get.EF2k <- function(f, k) {
  EF2k <- lapply(f[["EF2"]], function(X) X[, k])
  class(EF2k) <- "r1"
  return(EF2k)
}
get.KLk <- function(f, k) sapply(f[["KL"]], getElement, k)

get.n.fixed <- function(f)
  return(length(f[["fix.dim"]]))
get.fix.dim <- function(f, k = NULL) {
  if (is.null(f[["fix.dim"]]))
    return(NULL)
  if (is.null(k))
    return(f[["fix.dim"]])
  return(f[["fix.dim"]][[k]])
}
get.fix.idx <- function(f, k = NULL) {
  if (is.null(f[["fix.idx"]]))
    return(NULL)
  if (is.null(k))
    return(f[["fix.idx"]])
  return(f[["fix.idx"]][[k]])
}
get.fix.vals <- function(f, k = NULL) {
  if (is.null(f[["fix.vals"]]))
    return(NULL)
  if (is.null(k))
    return(f[["fix.vals"]])
  return(f[["fix.vals"]][[k]])
}

is.next.fixed <- function(f) {
  return(!is.null(get.next.fix.dim(f)))
}
get.next.fix.dim <- function(f) {
  if (is.null(f[["fix.dim"]]))
    return(NULL)
  next.k <- get.n.factors(f) + 1
  return(f[["fix.dim"]][[next.k]])
}
get.next.fix.idx <- function(f) {
  if (is.null(f[["fix.idx"]]))
    return(NULL)
  next.k <- get.n.factors(f) + 1
  return(f[["fix.idx"]][[next.k]])
}
get.next.fix.vals <- function(f) {
  if (is.null(f[["fix.vals"]]))
    return(NULL)
  next.k <- get.n.factors(f) + 1
  return(f[["fix.vals"]][[next.k]])
}
get.next.nonneg.dims <- function(f) {
  if (is.null(f[["nonneg.dims"]]))
    return(NULL)
  next.k <- get.n.factors(f) + 1
  return(f[["nonneg.dims"]][[next.k]])
}

get.next.unfixed.idx <- function(f) {
  next.k <- get.n.factors(f) + 1
  return(get.unfixed.idx(f, next.k))
}
get.unfixed.idx <- function(f, k) {
  fix.dim <- get.fix.dim(f, k)
  fix.idx <- get.fix.idx(f, k)
  return(setdiff(1:(get.dims(f)[[fix.dim]]), fix.idx))
}
all.fixed <- function(f, n) {
  fix.dim <- as.integer(get.fix.dim(f))
  idx.subset <- get.idx.subset(f)
  return(identical(fix.dim, n) && length(idx.subset) == 0)
}

add.subset.data <- function(factor, flash, fix.dim, idx.subset) {
  factor[["subset.data"]] <- get.subset.data(flash, fix.dim, idx.subset)
  factor[["idx.subset"]]  <- NULL
  return(factor)
}
get.subset.data <- function(f, fix.dim, idx.subset) {
  if (length(idx.subset) < 1)
    return(NULL)
  subset.data <- list(idx.subset = idx.subset)
  subset.data$R.subset  <- fullrank.subset(f[["R"]], fix.dim, idx.subset)
  subset.data$Y.subset  <- fullrank.subset(f[["Y"]], fix.dim, idx.subset)
  subset.data$Z.subset  <- fullrank.subset(f[["Z"]], fix.dim, idx.subset)
  subset.data$EF.subset <- lowrank.subset(f[["EF"]], fix.dim, idx.subset)
  return(subset.data)
}
get.idx.subset   <- function(f) {
  if (!is.null(f[["subset.data"]]))
    return(f[["subset.data"]][["idx.subset"]])
  return(f[["idx.subset"]])
}

set.R <- function(f, R) {
  f[["R"]] <- R
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
set.KL <- function(f, KL, n = NULL) {
  if (is.null(n)) {
    f[["KL"]] <- KL
  } else {
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
set.est.tau <- function(f, est.tau) {
  f[["est.tau"]] <- est.tau
  return(f)
}
set.R2 <- function(f, R2) {
  f[["R2"]] <- R2
  return(f)
}
set.obj <- function(f, obj) {
  f[["obj"]] <- obj
  return(f)
}
set.delta.R2 <- function(f, delta.R2) {
  f[["delta.R2"]] <- delta.R2
  return(f)
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

add.factor.to.EF <- function(f, new.EF) {
  if (is.null(f[["EF"]])) {
    f[["EF"]] <- as.lowrank(new.EF)
  } else {
    f[["EF"]] <- mapply(cbind, f[["EF"]], new.EF)
  }
  return(f)
}
add.factor.to.EF2 <- function(f, new.EF2) {
  if (is.null(f[["EF2"]])) {
    f[["EF2"]] <- as.lowrank(new.EF2)
  } else {
    f[["EF2"]] <- mapply(cbind, f[["EF2"]], new.EF2)
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
add.factor.to.g <- function(f, new.g) {
  f[["g"]] <- c(f[["g"]], list(new.g))
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
greedy.failed <- function(f) {
  return(identical(f[["greedy.fail"]], TRUE))
}

to.flashr <- function(f) {
  flash <- list()
  flash$EL     <- f$EF[[1]]
  flash$EF     <- f$EF[[2]]
  flash$EL2    <- f$EF2[[1]]
  flash$EF2    <- f$EF2[[2]]
  flash$fixl   <- matrix(FALSE, nrow = nrow(flash$EL), ncol = ncol(flash$EL))
  flash$fixf   <- matrix(FALSE, nrow = nrow(flash$EF), ncol = ncol(flash$EF))
  flash$KL_l   <- as.list(f$KL[[1]])
  flash$KL_f   <- as.list(f$KL[[2]])
  flash$tau    <- f$est.tau
  class(flash) <- "flash_fit"

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

  return(flash)
}
