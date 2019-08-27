#' @export
fitted.flash <- function(x) {
  if (x$n.factors == 0) {
    stop("Flash object does not have any factors.")
  }
  if (length(x$loadings.pm) == 2) {
    if (length(x$loadings.scale) == 1) {
      return(x$loadings.scale * x$loadings.pm[[1]] %*% t(x$loadings.pm[[2]]))
    } else {
      return(x$loadings.pm[[1]] %*% diag(x$loadings.scale) %*% t(x$loadings.pm[[2]]))
    }
  } else {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }
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
