#' Initialize flash factors at specified values
#'
#' Initializes factors at values specified by \code{init}. This
#'   has two primary uses: 1. One can initialize multiple factors at once using
#'   an SVD-like function and then optimize them via \code{flash.backfit}.
#'   Sometimes this results in a better fit than adding factors one at a time
#'   via \code{flash.add.greedy}. 2. One can initialize factors and then fix
#'   the loadings via \code{\link{flash.fix.loadings}} to incorporate "known"
#'   factors into a flash object. See below for examples of both use cases.
#'
#' @inheritParams flash
#'
#' @param flash A \code{flash} or \code{flash.fit} object to which factors are
#'   to be added.
#'
#' @param init An SVD-like object (specifically, a list containing fields
#'   \code{u}, \code{d}, and \code{v}), a flash fit, or a list of matrices
#'   specifying the values at which factors are to be initialized (for a data
#'   matrix of size \eqn{n \times p}, this should be a list of length two,
#'   with the first element a matrix of size \eqn{n \times k} and the second
#'   a matrix of size \eqn{p \times k}. If a flash fit is supplied, then it
#'   will be used to initialize both the first and second moments of
#'   posteriors on loadings and factors. Otherwise, the supplied values will
#'   be used to initialize posterior means, with posterior second moments
#'   initialized as the squared values of the first moments. Missing entries
#'   are not allowed.
#'
#' @examples
#' # Initialize several factors at once and backfit.
#' fl <- flash.init(gtex) %>%
#'   flash.init.factors(init = svd(gtex, nu = 5, nv = 5)) %>%
#'   flash.backfit()
#'
#' # Add a fixed factor with row loadings identically equal to one. This can be
#' #   interpreted as a "mean" factor that accounts for different row-wise means.
#' ones <- matrix(1, nrow = nrow(gtex), ncol = 1)
#' # Initialize the column loadings at the least squares solution.
#' ls.soln <- t(solve(crossprod(ones), crossprod(ones, gtex)))
#' fl <- flash.init(gtex) %>%
#'   flash.init.factors(init = list(ones, ls.soln)) %>%
#'   flash.fix.loadings(kset = 1, mode = 1L) %>%
#'   flash.backfit() %>%
#'   flash.add.greedy(Kmax = 5L)
#'
#' @importFrom ebnm ebnm_point_normal
#'
#' @export
#'
flash.init.factors <- function(flash,
                               init,
                               ebnm.fn = ebnm::ebnm_point_normal) {
  flash <- get.fit(flash)

  if (inherits(init, "flash")) {
    init <- get.fit(init)
  }

  EF2 <- NULL
  if (is.udv(init)) {
    EF <- udv.to.lowrank(init)
  } else if (inherits(init, "flash.fit")) {
    EF <- get.EF(init)
    EF2 <- get.EF2(init)
  } else if (is.list(init) && all(sapply(init, is.matrix))) {
    EF <- init
    class(EF) <- list("lowrank", "list")
  } else {
    stop("init must be an SVD-like object, a flash fit, or a list of matrices.")
  }

  if (is.null(EF2)) {
    EF2 <- lowrank.square(EF)
  }

  dims.must.match(EF, get.Y(flash))
  dims.must.match(EF2, get.Y(flash))

  if (anyNA(unlist(EF))) {
    stop("The initialization may not have missing data.")
  }

  ebnm.fn <- handle.ebnm.fn(ebnm.fn, get.dim(flash))

  flash <- set.EF(flash, lowranks.combine(get.EF(flash), EF))
  flash <- set.EF2(flash, lowranks.combine(get.EF2(flash), EF2))

  if (uses.R(flash)) {
    R <- get.Y(flash) - lowrank.expand(get.EF(flash))
    flash <- set.R(flash, get.nonmissing(flash) * R)
  }

  flash <- init.tau(flash)
  flash <- set.obj(flash, calc.obj(flash))

  K <- ncol(EF[[1]])

  # Initialize KL at zero and g at NULL.
  for (n in 1:get.dim(flash)) {
    flash <- set.KL(flash, c(get.KL(flash, n), rep(0, K)), n)
  }
  EF.g <- rep(list(rep(list(NULL), get.dim(flash))), K)
  flash <- set.g(flash, c(get.g(flash), EF.g))

  ebnm.fn <- c(get.ebnm.fn(flash), rep(list(ebnm.fn), length.out = K))
  flash <- set.ebnm.fn(flash, ebnm.fn)

  # Initialize is.valid and is.zero.
  flash <- set.is.valid(flash, c(is.valid(flash), rep(FALSE, K)))
  flash <- set.is.zero(flash, c(is.zero(flash), rep(FALSE, K)))

  flash <- wrapup.flash(flash, output.lvl = 3L)

  return(flash)
}
