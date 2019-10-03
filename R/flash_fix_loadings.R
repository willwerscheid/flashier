# TODO: document and export
# TODO: really should be able to fix more than one mode, no?
# allow is.fixed to be matrix or vector (when same for each), but must be logical
flash.fix.loadings <- function(flash, kset, mode, is.fixed = TRUE) {
  # TODO: argument checking
  is.fixed <- matrix(is.fixed, nrow = get.dims(flash)[mode], ncol = length(kset))

  fix.dim <- get.fix.dim(flash)
  fix.idx <- get.fix.idx(flash)

  # TODO: check if already set; if mode different throw error
  for (i in 1:length(kset)) {
    k <- kset[i]

    if ((length(fix.dim) >= k)
        && !is.null(fix.dim[[k]])
        && (fix.dim[[k]] != mode)) {
      stop(paste("Loadings can only be fixed along a single mode for any",
                 "given factor."))
    }

    fix.dim[[k]] <- mode
    fix.idx[[k]] <- which(is.fixed[, i])
  }

  flash <- set.fix.dim(flash, fix.dim)
  flash <- set.fix.idx(flash, fix.idx)

  return(flash)
}

# TODO: document and export
flash.unfix.loadings <- function(flash, kset) {
  fix.dim <- get.fix.dim(flash)
  fix.idx <- get.fix.idx(flash)

  for (k in kset) {
    fix.dim[[k]] <- NULL
    fix.idx[[k]] <- NULL
  }

  flash <- set.fix.dim(flash, fix.dim)
  flash <- set.fix.idx(flash, fix.idx)

  return(flash)
}
