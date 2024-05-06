#' Set timeout
#'
#' Used in a \code{\link{flash}} pipeline to set a maximum amount of fitting
#'   time. Note that timeout conditions are only checked during greedy fits
#'   and backfits, so that the total amount of fitting time can exceed the time
#'   set by \code{flash_set_timeout} (especially if, for example, there is a
#'   nullcheck involving many factor/loading pairs). Also note that timeout
#'   conditions must be cleared using function \code{\link{flash_clear_timeout}}
#'   before any re-fitting is attempted.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param tim A numeric value giving the maximum amount of fitting time, with
#'   the units of time specified by parameter \code{units}.
#'
#' @param units The units of time according to which parameter \code{tim} is to
#'   be interpreted.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with the
#'   timeout settings reflected in updates to the "internal" \code{flash_fit}
#'   object. These settings will persist across all subsequent calls to
#'   \code{flash_xxx} functions until they are modified either by
#'   \code{\link{flash_clear_timeout}} or by another call to
#'   \code{flash_set_timeout}.
#'
#' @examples
#' fl <- flash_init(gtex) |>
#'   flash_set_timeout(1, "secs") |>
#'   flash_greedy(Kmax = 30) |>
#'   flash_backfit() |>
#'   flash_nullcheck() |>
#'   flash_clear_timeout() # Always clear timeout at the end of a pipeline.
#'
#' @export
#'
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

#' Set timeout
#'
#' Used in a \code{\link{flash}} pipeline to clear timeout conditions set using
#'   \code{\link{flash_set_timeout}}.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   timeout settings cleared.
#'
#' @export
#'
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
