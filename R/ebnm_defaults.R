# Output required for usual flash updates.
default.ebnm.output <- c("posterior_mean",
                         "posterior_second_moment",
                         "fitted_g",
                         "log_likelihood")

# Since ebnm will typically be called many times, it suffices to do a small
#   number of optimization iterations each time.
mixsqp.defaults <- list(maxiter.sqp = 10)
nlm.defaults <- list(iterlim = 10)

# Since maxiter.sqp is intentionally set to be small, ignore the ashr warning
#   about reaching the maximum number of iterations.
# Also ignore the ebnm warning about setting mode and scale when g is fixed.
#   This pops up when calculating LFSR and posterior samplers using
#   nonzero.mode prior families.
ignored.warnings <- c("Optimization failed to converge. Results",
                      "mode and scale parameters are ignored")

#' @importFrom ebnm ebnm
ebnm.nowarn <- function(...) {
  withCallingHandlers(res <- ebnm(...),
                      warning = function(w) {
                        if (any(startsWith(conditionMessage(w), ignored.warnings)))
                          invokeRestart("muffleWarning")
                      })
  return(res)
}
