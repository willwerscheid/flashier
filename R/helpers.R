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

get.EF <- function(f) {
  EF <- f[["EF"]]
  if (is.null(EF[[1]]))
    return(NULL)
  return(EF)
}
get.EF2 <- function(f) {
  EF2 <- f[["EF2"]]
  if (is.null(EF2[[1]]))
    return(NULL)
  return(EF2)
}

uses.R <- function(f) !is.null(f[["R"]])

is.zero <- function(f, k = NULL) {
  if (is.null(k))
    return(f[["is.zero"]])
  return(f[["is.zero"]][k])
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
set.to.zero <- function(f, k = NULL) {
  if (is.null(k)) {
    f[["is.zero"]] <- TRUE
  } else {
    f[["is.zero"]][k] <- TRUE
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
  return(flash)
}
