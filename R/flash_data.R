set.flash.data <- function(data, S = NULL, S.dim = NULL) {
  must.be.supported.data.type(data, allow.null = FALSE)
  must.be.supported.data.type(S, allow.vector = TRUE)
  must.be.integer(S.dim, lower = 0, upper = length(dim(data)))

  flash.data <- list()
  flash.data$Y <- data
  if (anyNA(data)) {
    flash.data$Y[is.na(data)] <- 0
    flash.data$Z <- 1L * !is.na(data)
  } else {
    flash.data$Z <- 1
  }

  if (is.vector(S) && is.null(S.dim)) {
    if (length(S) == 1) {
      S.dim <- 0
    } else {
      S.dim <- which(length(S) == dim(data))
      if (length(S.dim) == 0)
        stop("S was interpreted as a vector, but couldn't be aligned ",
             "with the data.")
      if (length(S.dim) > 1)
        stop("S could not be unambiguously interpreted. Set data using ",
             "set.flash.data with S.dim specified.")
    }
  } else {
    dims.must.match(data, S, S.dim)
  }

  if (!is.null(S)) {
    tau <- 1 / S^2
  } else {
    tau <- NULL
  }

  flash.data$given.tau     <- tau
  flash.data$given.tau.dim <- S.dim

  class(flash.data) <- "flash.data"

  return(flash.data)
}
