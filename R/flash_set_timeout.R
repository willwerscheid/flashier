flash_set_timeout <- function(flash,
                              tim,
                              units = c("hours", "secs", "mins", "days", "weeks")) {
  fit <- get.fit(flash)

  must.be.numeric(tim, allow.null = FALSE)
  units <- match.arg(units)
  timeout <- Sys.time() + as.difftime(tim, units = units)

  fit <- clear.timeout.reached.flag(fit)
  fit <- set.timeout(fit, timeout)

  flash <- set.fit(flash, fit)

  return(flash)
}

flash_clear_timeout <- function(flash) {
  fit <- get.fit(flash)

  fit <- clear.timeout.reached.flag(fit)
  fit <- set.timeout(fit, NULL)

  flash <- set.fit(flash, fit)

  return(flash)
}

is.timed.out <- function(f) {
  return(!is.null(get.timeout(f)) && get.timeout(f) < Sys.time())
}
