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
#' @return The \code{\link{flash}} object from argument \code{flash}, with
#'   the data modified as specified by \code{newdata}. Residual variances and
#'   ELBO are also updated.
#'
#' @export
#'
flash_update_data <- function(flash, newdata) {
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
  flash <- set.Y2(flash, NULL) # Will be recomputed by init.tau().
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
