set.default.tol <- function(flash) {
  return(sqrt(.Machine$double.eps) * prod(get.dims(flash)))
}

get.conv.crit <- function(update.info) {
  return(update.info[[length(update.info)]])
}
