calc.update.info <- function(new, old, conv.crit.fn, verbose.fns, k = NULL) {
  if (length(verbose.fns) == 0) {
    all.fns <- list(conv.crit.fn)
  } else {
    all.fns <- verbose.fns
    conv.crit.idx <- which(sapply(verbose.fns, identical, conv.crit.fn))
    if (length(conv.crit.idx) == 0)
      all.fns <- c(all.fns, conv.crit.fn)
  }

  update.info <- lapply(all.fns, do.call, list(new, old, k))

  # Put a copy of the convergence criterion at the end of the list so that
  #   it's easy to find.
  if ((length(verbose.fns) > 0) && (length(conv.crit.idx) > 0))
    update.info <- c(update.info, update.info[conv.crit.idx])

  return(update.info)
}

get.verbose.info <- function(update.info) {
  return(update.info[-length(update.info)])
}

get.new.obj <- function(new, old, k) {
  if (!is.obj.valid(new))
    return(NA)
  return(get.obj(new))
}

display.obj <- function(new, old, k) {
  obj <- get.new.obj(new, old, k)
  if (is.na(obj))
    return("NA")
  return(formatC(obj, format = "f", digits = 2))
}

calc.obj.diff <- function(new, old, k) {
  if (!is.obj.valid(old) || !is.obj.valid(new))
    return(Inf)
  return(get.obj(new) - get.obj(old))
}

calc.max.abs.chg.EF <- function(new, old, k, n = NULL) {
  if (!is.null(k))
    return(calc.max.abs.chg(get.EF.k(new, k, n), get.EF.k(old, k, n)))
  return(calc.max.abs.chg(get.EF(new, n), get.EF(old, n)))
}

calc.max.chg.EF <- function(new, old, k, n = NULL) {
  if (!is.null(k))
    return(calc.max.chg(get.EF.k(new, k, n), get.EF.k(old, k, n)))
  return(calc.max.chg(get.EF(new, n), get.EF(old, n)))
}

which.max.chg.EF <- function(new, old, k, n = NULL) {
  if (!is.null(k))
    return(which.max.chg(get.EF.k(new, k, n), get.EF.k(old, k, n)))
  return(which.max.chg(get.EF(new, n), get.EF(old, n)))
}

get.sparsity <- function(new, old, k, n) {
  if (!is.null(k)) {
    g <- get.g.k(new, k, n)
  } else {
    g <- get.g(new, n)
  }
  return(g$pi[1])
}

get.exclusion.count <- function(new, old, k, n) {
  if (!is.null(k)) {
    return(length(get.exclusions(new)[[k]][[n]]))
  } else {
    return(length(get.exclusions(new)[[n]]))
  }
}

calc.max.abs.chg <- function(new, old) {
  new <- l2.normalize.and.stack(new)
  old <- l2.normalize.and.stack(old)
  return(max(abs(new - old)))
}

calc.max.chg <- function(new, old) {
  new <- l2.normalize.and.stack(new)
  old <- l2.normalize.and.stack(old)
  max.increase <- max(new - old)
  max.decrease <- max(old - new)
  if (max.increase > max.decrease)
    return(max.increase)
  return(-max.decrease)
}

which.max.chg <- function(new, old) {
  new <- l2.normalize.and.stack(new)
  old <- l2.normalize.and.stack(old)
  return(which.max(apply(abs(new - old), 1, max)))
}

l2.normalize.and.stack <- function(x) {
  if (is.list(x)) {
    norm.x <- lapply(x, l2.normalize)
    return(do.call(rbind, norm.x))
  }
  return(l2.normalize(x))
}

l2.normalize <- function(x) {
  if (is.matrix(x)) {
    norm <- sqrt(colSums(x^2))
  } else {
    norm <- sqrt(sum(x^2))
  }
  norm[norm == 0] <- 1
  if (is.matrix(x))
    return(x / matrix(norm, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))
  return(matrix(x / norm, ncol = 1))
}
