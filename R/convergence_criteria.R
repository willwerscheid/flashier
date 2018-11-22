# TODO: more functions for verbose output: sparsity, max chg for single dim

calc.update.info <- function(new, old, conv.crit.fn, verbose.fns) {
  if (length(verbose.fns) == 0) {
    all.fns <- list(conv.crit.fn)
  } else {
    all.fns <- verbose.fns
    conv.crit.idx <- which(sapply(verbose.fns, identical, conv.crit.fn))
    if (length(conv.crit.idx) == 0)
      all.fns <- c(all.fns, conv.crit.fn)
  }

  update.info <- sapply(all.fns, do.call, list(new, old))

  if (length(verbose.fns) > 0 && length(conv.crit.idx) > 0)
    update.info <- c(update.info, update.info[conv.crit.idx])

  return(update.info)
}

get.conv.crit <- function(update.info) {
  return(update.info[length(update.info)])
}

calc.obj.diff <- function(new, old) {
  if (!is.obj.valid(old) || !is.obj.valid(new))
    return(Inf)
  return(get.obj(new) - get.obj(old))
}

calc.max.chg <- function(new, old) {
  if (is.list(new)) {
    new <- lapply(new, l2.normalize)
    old <- lapply(old, l2.normalize)
  } else {
    new <- l2.normalize(new)
    old <- l2.normalize(old)
  }
  return(max(abs(unlist(new) - unlist(old))))
}

l2.normalize <- function(x) {
  norm <- sum(x^2)
  if (norm == 0)
    return(x)
  return(x / sum(x^2))
}
