#' @export
fitted.flash <- function(x) {
  if (x$n.factors == 0) {
    stop("Flash object does not have any factors.")
  }
  if (!is.null(x$factors.pm)) {
    return(x$loadings.pm %*% t(x$factors.pm))
  } else {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }
}

#' @export
fitted.flash.fit <- function(x) {
  if (uses.R(x)) {
    R <- get.R(x)
  } else {
    R <- get.Y(x) - lowrank.expand(get.EF(x))
  }

  if (any.missing(x)) {
    R[get.nonmissing(x) == 0] <- NA
  }

  return(R)
}

#' @export
ldf <- function(x, ...) {
  UseMethod("ldf", x)
}

#' @export
ldf.flash <- function(x, type = "f") {
  return(ldf.flash.fit(get.fit(x), type = type))
}

#' @export
ldf.flash.fit <- function(x, type = "f") {
  type <- tolower(type)
  if (type == "2") {
    type <- "f"
  } else if (type == "1") {
    type <- "o"
  } else if (type == "m") {
    type <- "i"
  }

  if (get.n.factors(x) == 0) {
    stop("Flash fit does not have any factors.")
  }
  if (get.dim(x) > 2) {
    stop("S3 method \"ldf\" not available for tensors.")
  }

  ldf <- calc.normalized.loadings(x, type = type)

  ret <- list()
  ret$L <- ldf$normalized.loadings[[1]]
  ret$D <- ldf$scale.constants
  ret$F <- ldf$normalized.loadings[[2]]

  return(ret)
}

calc.normalized.loadings <- function(flash, for.pve = FALSE, type = "f") {
  ret <- list()

  if (for.pve) {
    loadings <- get.EF2(flash)
    norms <- lapply(loadings, colSums)
  } else {
    loadings <- get.EF(flash)
    if (type == "f") {
      norms <- lapply(loadings, function(x) {sqrt(colSums(x^2))})
    } else if (type == "o") {
      norms <- lapply(loadings, function(x) {colSums(abs(x))})
    } else if (type == "i") {
      norms <- lapply(loadings, function(x) {apply(abs(x), 2, max)})
    } else {
      stop("Norm type not recognized.")
    }
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

#' @export
print.flash = function(x, ...) {
  if (x$n.factors == 0) {
    cat("Flash object with zero factors.\n")
  } else if (x$n.factors == 1) {
    cat("Flash object with one factor.\n")
    cat(sprintf("  Proportion of variance explained: %0.3f\n", x$pve))
  } else {
    cat(sprintf("Flash object with %d factors.\n", x$n.factors))
    if (all(x$pve < 0.001)) {
      cat("  All factors have PVE < 0.001.")
    } else {
      cat("  Proportion of variance explained")
      if (any(x$pve < 0.001)) {
        cat("*")
      }
      cat(":\n")
      for (k in 1:length(x$pve)) {
        if (x$pve[k] > 0.001) {
          cat(sprintf("    Factor %d: %0.3f\n", k, x$pve[k]))
        }
      }
      if (any(x$pve < 0.001)) {
        cat(paste("    *Factors with PVE < 0.001 are omitted from this",
                  "summary.\n"))
      }
    }
  }

  if (!is.na(x$elbo)) {
    cat(sprintf("  Variational lower bound: %0.3f\n", x$elbo))
  }
}
