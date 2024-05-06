#' Add "intercept" to a flash object
#'
#' Adds an all-ones vector as a fixed set of loadings (if \code{rowwise = TRUE})
#'   or fixed factor (if \code{rowwise = FALSE}). Assuming (without loss of
#'   generality) that the fixed factor/loadings is indexed as \eqn{k = 1},
#'   a fixed set of loadings gives:
#'   \deqn{\mathbf{y}_{i \cdot} \approx \mathbf{f}_1 + \sum_{k = 2}^K \ell_{i k}
#'   \mathbf{f}_k,}
#'   so that the (estimated) factor \eqn{\mathbf{f}_1 \in \mathbf{R}^p} is shared
#'   by all row-wise observations \eqn{\mathbf{y}_{i \cdot} \in \mathbf{R}^p}.
#'   A fixed factor gives:
#'   \deqn{\mathbf{y}_{\cdot j} \approx \boldsymbol{\ell}_1 + \sum_{k = 2}^K f_{j k}
#'   \boldsymbol{\ell}_k,}
#'   so that the (estimated) set of loadings \eqn{\ell_1 \in \mathbf{R}^n} is
#'   shared by all column-wise observations \eqn{y_{\cdot j} \in \mathbf{R}^n}.
#'
#'   The estimated factor (if \code{rowwise = TRUE}) or set of loadings
#'   (if \code{rowwise = FALSE}) is initialized at the column-
#'   or row-wise means of the data (or, if factor/loadings pairs have previously
#'   been added, at the column- or row-wise means of the matrix of residuals)
#'   and then backfit via function \code{\link{flash_backfit}}.
#'
#' @param flash A \code{flash} or \code{flash_fit} object to which an "intercept"
#'   is to be added.
#'
#' @param rowwise Should the all-ones vector be added as a fixed set of loadings
#'   ("row-wise") or a fixed factor ("column-wise")? See above for details.
#'
#' @param ebnm_fn As with other factor/loadings pairs, a prior is put on the
#'   estimated factor (if \code{rowwise = TRUE}) or set of loadings (if
#'   \code{rowwise = FALSE}). Parameter \code{ebnm_fn} specifies the function
#'   used to estimate that prior; see \code{\link{flash}} for details.
#'
#' @examples
#' # The following are equivalent:
#' init <- list(matrix(rowMeans(gtex), ncol = 1),
#'              matrix(1, nrow = ncol(gtex)))
#' fl <- flash_init(gtex) |>
#'   flash_factors_init(init) |>
#'   flash_factors_fix(kset = 1, which_dim = "factors") |>
#'   flash_backfit(kset = 1)
#'
#' fl <- flash_init(gtex) |>
#'   flash_add_intercept(rowwise = FALSE)
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with an
#'   "intercept" added.
#'
#' @importFrom stats residuals
#'
#' @export
#'
flash_add_intercept <- function(flash,
                                rowwise = TRUE,
                                ebnm_fn = ebnm_point_normal) {
  fit <- flash_fit(flash)

  if (rowwise) {
    fixed_dim <- "loadings"
    ones <- matrix(1, nrow = get.dims(fit)[1])
    if (any_missing(fit)) {
      init.F <- colSums(residuals(fit), na.rm = TRUE) / colSums(get.nonmissing(fit))
    } else {
      init.F <- colMeans(residuals(fit))
    }
    init <- list(ones, matrix(init.F, ncol = 1))
  } else {
    fixed_dim <- "factors"
    ones <- matrix(1, nrow = get.dims(fit)[2])
    if (any_missing(fit)) {
      init.L <- rowSums(residuals(fit), na.rm = TRUE) / rowSums(get.nonmissing(fit))
    } else {
      init.L <- rowMeans(residuals(fit))
    }
    init <- list(matrix(init.L, ncol = 1), ones)
  }

  flash <- flash |>
    flash_factors_init(init, ebnm_fn = ebnm_fn) |>
    flash_factors_fix(kset = get.n.factors(fit) + 1, which_dim = fixed_dim) |>
    flash_backfit(kset = get.n.factors(fit) + 1)

  return(flash)
}
