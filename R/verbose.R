# Level 1 announcements (main loop) -------------------------------------------

announce.lvl1 <- function(verbose.lvl) {
  return(verbose.lvl > 0)
}

announce.flash.init <- function(verbose.lvl) {
  if (announce.lvl1(verbose.lvl))
    message("Initializing flash object...")
}

announce.greedy <- function(verbose.lvl, k) {
  if (announce.lvl1(verbose.lvl))
    message("Adding factor ", k, " to flash object...")
}

report.greedy.result <- function(verbose.lvl, failure) {
  if (announce.lvl1(verbose.lvl) && failure) {
    message("Failed to add new factor.")
  }
}

announce.backfit <- function(verbose.lvl, n.factors) {
  if (announce.lvl1(verbose.lvl))
    message("Backfitting ", n.factors, " factors...")
}

announce.nullchk <- function(verbose.lvl, n.factors) {
  if (announce.lvl1(verbose.lvl))
    message("Nullchecking ", n.factors, " factors...")
}

announce.wrapup <- function(verbose.lvl) {
  if (announce.lvl1(verbose.lvl))
    message("Wrapping up...")
}

report.completion <- function(verbose.lvl) {
  if (announce.lvl1(verbose.lvl))
    message("Done.")
}

# Level 2 announcements -------------------------------------------------------

announce.lvl2 <- function(verbose.lvl) {
  return(verbose.lvl > 1)
}

announce.factor.init <- function(verbose.lvl) {
  if (announce.lvl2(verbose.lvl))
    message("  Initializing factor...")
}

announce.factor.opt <- function(verbose.lvl) {
  if (announce.lvl2(verbose.lvl))
    message("  Optimizing factor...")
}

# Optimization details (level 3) ----------------------------------------------
