# Check that an initial value of g has a fighting chance of fitting the data.
#
#' @importFrom stats pnorm
#'
warmstart.sanity.check = function(g, x, s) {
  # Find rough limits for the region where g can have significantly positive
  #   density (pi can be ignored because it will be re-estimated).
  if (inherits(g,"unimix")) {
    upper.grange = max(g$b)
    lower.grange = min(g$a)
  } else if (inherits(g, "normalmix")) {
    # In the normal and halfnormal cases, use an anti-conservative range. It's
    #   better to re-estimate the grid than to use a bad one.
    upper.grange = 2 * max(g$sd)
    lower.grange = -upper.grange
  } else if (inherits(g, "tnormalmix")) {
    upper.grange = 2 * max(g$sd[is.infinite(g$b)])
    lower.grange = -2 * max(g$sd[is.infinite(g$a)])
  } else {
    # Skip the sanity check if g is not from a recognized ashr class.
    return(TRUE)
  }

  # Find the most outlying data points and calculate p-values conditional
  #   on the true values lying at the limits of g's range.
  worst.lower = min((x - lower.grange) / s)
  lower.p = pnorm(worst.lower)
  worst.upper = max((x - upper.grange) / s)
  upper.p = 1 - pnorm(worst.upper)

  # In the +uniform case, the worst lower point is ignored because we can't
  #   do anything about it.
  if (length(g$a) > 1 && min(g$a) == max(g$a)) {
    lower.p = 1
  }
  # And similarly for the -uniform case.
  if (length(g$b) > 1 && min(g$b) == max(g$b)) {
    upper.p = 1
  }

  worst.p = min(lower.p, upper.p)

  # The case where g is a single (null) component and the data comes from
  #   the null model should pass with 95% probability.
  n = length(x)
  passes.sanity.check = (worst.p > 1 - 0.95^(1/n))

  return(passes.sanity.check)
}
