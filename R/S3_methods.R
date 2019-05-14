fitted.flash <- function(x) {
  if (get.dim(x$fit)) {
    x$fit$EF[[1]] %*% t(x$fit$EF[[2]])
  } else {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }
}
