#' Update data in a flash object
#'
#' Replaces the data in a flash object with a new set of observations. Estimates
#'   of residual variances and the ELBO are also updated.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param newdata The new observations. Can be a matrix, a sparse matrix of
#'   class \code{\link[Matrix]{Matrix}}, or a low-rank matrix representation.
#'
#' @param Y2_diff Optionally, users can supply the (summed) changes in the
#'   squared values of the data \eqn{y_{ij}^2}, which are needed to estimate the
#'   residual variance parameters \eqn{s_{ij}^2} for simple variance structures
#'   (i.e., when \code{var_type} is set to 0, 1, or 2). If calculating
#'   entries \eqn{y_{ij}^2} from scratch is expensive, supplying an argument
#'   to \code{Y2_diff} can greatly speed up data updates. If specified, the
#'   argument should be a scalar
#'   \eqn{\sum_{i, j} \left( y_{ij}^{2 \text{(new)}} - y_{ij}^{2 \text{(old)}} \right)}
#'   when \code{var_type = 0}; a vector of length \eqn{n} with entries
#'   \eqn{\sum_{j = 1}^p \left( y_{ij}^{2 \text{(new)}} - y_{ij}^{2 \text{(old)}} \right)}
#'   when \code{var_type = 1}; or a vector of length \eqn{p} with entries
#'   \eqn{\sum_{i = 1}^n \left( y_{ij}^{2 \text{(new)}} - y_{ij}^{2 \text{(old)}} \right)}
#'   when \code{var_type = 2}. The argument is ignored when any other variance
#'   structure is used.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   the data modified as specified by \code{newdata}. Residual variances and
#'   ELBO are also updated.
#'
#' @export
#'
flash_update_data <- function(flash, newdata, Y2_diff = NULL) {
  flash <- get.fit(flash)
  fl.dims <- get.dims(flash)

  # Returns "flash.data" object with fields Y and Z:
  data <- set.flash.data(newdata)
  data.dims <- get.dims(data)

  if (!identical(data.dims, fl.dims)) {
    stop("Dimensions of newdata must match dimensions of existing flash data.")
  }

  # Update residuals:
  if (is.var.type.zero(flash) && !is.tau.simple(flash)) {
    old.data <- get.Y(flash, require.fullrank = TRUE)
    flash <- set.R(flash, get.R(flash) + get.Y(data) - old.data)
  }

  # Update data (Y and Z):
  Z.changed <- identical(get.nonmissing(data), get.nonmissing(flash))
  flash <- set.Y(flash, get.Y(data))
  if (is.null(Y2_diff) || is.null(get.Y2(flash))) {
    flash <- set.Y2(flash, NULL) # Will be recomputed by init.tau().
  } else {
    if (length(Y2_diff) != length(get.Y2(flash))) {
      stop("Expected argument Y2_diff to have length ", length(get.Y2(flash)), ".")
    }
    flash <- set.Y2(flash, get.Y2(flash) + Y2_diff)
  }
  flash <- set.nonmissing(flash, get.nonmissing(data))

  # Update precomputed quantities:
  if (Z.changed) {
    if (is.tau.simple(flash)) {
      flash$n.nonmissing <- init.n.nonmissing(flash, get.R2.n(flash))
    } else if (is.var.type.kronecker(flash)) {
      flash$kron.nonmissing <- init.kron.nonmissing(flash)
    }
  }

  # Update residual variances and ELBO.
  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))

  flash <- wrapup.flash(flash, output.lvl = 3L)
  flash <- flash_set_verbose(flash, verbose = 1L)

  return(flash)
}
