optimize.it <- function(x,
                        update.fn,
                        update.args = NULL,
                        obj.fn,
                        obj.args = NULL,
                        tol,
                        maxiter = 100) {
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

calc.max.chg.r1 <- function(old, new) {
  old <- lapply(old, l2.normalize)
  new <- lapply(new, l2.normalize)
  return(max(abs(unlist(new) - unlist(old))))
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
