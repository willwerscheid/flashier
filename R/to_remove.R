add.next.factor <- function(flash,
                            greedy.tol = 1e-2,
                            greedy.maxiter = 200,
                            init.tol = 1e-2,
                            init.maxiter = 100,
                            init.fn = NULL) {
  factor <- init.factor(flash, init.fn, init.tol, init.maxiter, verbose = FALSE)
  if (is.fixed(factor)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- add.greedy(factor, flash, greedy.tol, greedy.maxiter)
  }
  return(flash)
}

add.greedy <- function(factor, flash, tol, maxiter) {
  factor <- optimize.it(factor,
                        update.fn = update.factor,
                        update.args = list(flash = flash),
                        obj.fn = calc.obj.diff,
                        tol = tol,
                        maxiter = maxiter)
  if (get.obj(factor) > get.obj(flash) || !is.obj.valid(flash, factor)) {
    flash <- add.new.factor.to.flash(factor, flash)
  } else {
    flash <- set.greedy.fail.flag(flash)
  }
  # TODO: what if objective decreases?

  return(flash)
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
    obj     <- do.call(obj.fn, c(list(x, old.x), obj.args))
    # TODO: verbose output here
    # TODO: handle obj decreases
    # message(obj)
  }

  return(x)
}

backfit <- function(flash, kset, shuffle.kset = FALSE, tol = 1e-2) {
  flash <- optimize.it(flash,
                       update.fn = backfit.once,
                       update.args = list(kset = kset,
                                          shuffle.kset = shuffle.kset),
                       obj.fn = calc.obj.diff,
                       tol = tol)
  return(flash)
}

backfit.once <- function(flash, kset, shuffle.kset = FALSE) {
  if (shuffle.kset)
    kset <- sample(kset)

  for (k in kset)
    flash <- update.kth.factor(flash, k)

  return(flash)
}
