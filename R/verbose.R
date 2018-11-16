# Level 1 announcements (main loop) -------------------------------------------

announce.flash.init <- function(verbose) {
  if (verbose)
    message("Initializing flash object...")
}

announce.greedy <- function(verbose) {
  if (verbose)
    message("Attempting to add factor to flash object...")
}

announce.backfit <- function(verbose) {
  if (verbose)
    message("Backfitting flash object...")
}

announce.nullchk <- function(verbose) {
  if (verbose)
    message("Performing nullcheck...")
}

# Level 2 announcements -------------------------------------------------------

announce.factor.init <- function(verbose) {
  if (verbose)
    message("  Initializing factor...")
}

announce.factor.opt <- function(verbose) {
  if (verbose)
    message("  Optimizing factor...")
}

# Optimization details (level 3) ----------------------------------------------
