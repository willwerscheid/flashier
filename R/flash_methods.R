#' @export
fitted.flash <- function(x) {
  f <- get.fit(x)
  if (get.dim(f)) {
    f$EF[[1]] %*% t(f$EF[[2]])
  } else {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }
}
