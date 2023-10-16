#' Calculate the difference in ELBO
#'
#' The default objective function used to determine convergence when fitting
#'   a \code{\link{flash}} object. Calculates the difference in the
#'   variational lower bound ("ELBO") from one iteration to the next.
#'
#' @details This function is an example of a function that may be passed to
#'   parameter \code{fn} in function \code{\link{flash_set_conv_crit}} to set
#'   the convergence criterion for a flash pipeline. See
#'   \code{\link{flash_set_conv_crit}} for details and examples.
#'
#' @inheritParams flash_verbose_elbo
#'
#' @return A scalar, which is compared against the tolerance parameter
#'   \code{tol} to determine whether a fit has converged.
#'
#' @seealso \code{\link{flash_conv_crit_max_chg}}
#'   \code{\link{flash_conv_crit_max_chg_L}},
#'   \code{\link{flash_conv_crit_max_chg_F}}
#'
#' @export
#'
flash_conv_crit_elbo_diff <- function(curr, prev, k) {
  return(calc.obj.diff(curr, prev, k))
}

#' Calculate the maximum absolute difference in scaled loadings and factors
#'
#' An alternative objective function that can be used to determine
#'   convergence when fitting a \code{\link{flash}} object. Calculates the
#'   maximum (absolute) change over all (posterior expected values for)
#'   loadings \eqn{\ell_{ik}} and factors \eqn{f_{jk}}. At each iteration, the
#'   loadings vectors \eqn{\ell_{\cdot 1}, \ldots, \ell_{\cdot K}} and factors
#'   \eqn{f_{\cdot 1}, \ldots, f_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_conv_crit_elbo_diff return
#'
#' @seealso \code{\link{flash_conv_crit_elbo_diff}},
#'   \code{\link{flash_conv_crit_max_chg_L}}
#'   \code{\link{flash_conv_crit_max_chg_F}}
#'
#' @export
#'
flash_conv_crit_max_chg <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = NULL))
}

#' Calculate the maximum absolute difference in scaled loadings
#'
#' An alternative objective function that can be used to determine
#'   convergence when fitting a \code{\link{flash}} object. Calculates the
#'   maximum (absolute) change over all (posterior expected values for)
#'   loadings \eqn{\ell_{ik}}. At each iteration, the loadings vectors
#'   \eqn{\ell_{\cdot 1}, \ldots, \ell_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_conv_crit_elbo_diff return
#'
#' @seealso \code{\link{flash_conv_crit_elbo_diff}},
#'   \code{\link{flash_conv_crit_max_chg}}
#'   \code{\link{flash_conv_crit_max_chg_F}}
#'
#' @export
#'
flash_conv_crit_max_chg_L <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = 1))
}

#' Calculate the maximum absolute difference in scaled factors
#'
#' An alternative objective function that can be used to determine
#'   convergence when fitting a \code{\link{flash}} object. Calculates the
#'   maximum (absolute) change over all (posterior expected values for)
#'   factors \eqn{f_{jk}}. At each iteration, the factors
#'   \eqn{f_{\cdot 1}, \ldots, f_{\cdot K}} are \eqn{L^2}-normalized.
#'
#' @inheritParams flash_verbose_elbo
#'
#' @inherit flash_conv_crit_elbo_diff return
#'
#' @seealso \code{\link{flash_conv_crit_elbo_diff}},
#'   \code{\link{flash_conv_crit_max_chg}}
#'   \code{\link{flash_conv_crit_max_chg_L}}
#'
#' @export
#'
flash_conv_crit_max_chg_F <- function(curr, prev, k) {
  return(calc.max.abs.chg.EF(curr, prev, k, n = 2))
}
