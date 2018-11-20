calc.update.info <- function(old, new, conv.crit.fn, verbose.fns) {
  if (length(verbose.fns) == 0) {
    all.fns <- list(conv.crit.fn)
  } else {
    all.fns <- verbose.fns
    conv.crit.idx <- which(sapply(verbose.fns, identical, conv.crit.fn))
    if (length(conv.crit.idx) == 0)
      all.fns <- c(all.fns, conv.crit.fn)
  }

  update.info <- sapply(all.fns, do.call, list(old, new))

  if (length(verbose.fns) > 0 && length(conv.crit.idx) > 0)
    update.info <- c(update.info, update.info[conv.crit.idx])

  return(update.info)
}

get.conv.crit <- function(update.info) {
  return(update.info[length(update.info)])
}

calc.max.chg.list <- function(old, new) {
  old <- lapply(old, l2.normalize)
  new <- lapply(new, l2.normalize)
  return(max(abs(unlist(new) - unlist(old))))
}

calc.max.chg <- function(old, new) {
  return(max(abs(unlist(l2.normalize(new)) - unlist(l2.normalize(old)))))
}

calc.max.EF.chg <- function(old, new) {
  return(calc.max.chg.list(get.EF(old), get.EF(new)))
}

simple.obj.diff <- function(old, new) {
  return(get.obj(new) - get.obj(old))
}

calc.obj.diff <- function(old, new) {
  if (!is.obj.valid(old) || !is.obj.valid(new))
    return(Inf)
  return(get.obj(new) - get.obj(old))
}

l2.normalize <- function(x) {
  norm <- sum(x^2)
  if (norm == 0)
    return(x)
  return(x / sum(x^2))
}

optimize.it <- function(x,
                        update.fn,
                        update.args = NULL,
                        obj.fn,
                        obj.args = NULL,
                        tol,
                        maxiter = 100,
                        verbose = NULL) {
  obj <- Inf
  iter <- 0
  while (obj > tol && iter < maxiter) {
    iter    <- iter + 1
    old.x   <- x
    x       <- do.call(update.fn, c(list(x), update.args))
    obj     <- do.call(obj.fn, c(list(old.x, x), obj.args))
    # TODO: verbose output here
    # TODO: handle obj decreases
    # message(obj)
  }

  return(x)
}

