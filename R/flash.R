#' Empirical Bayes matrix/tensor factorization
#'
#' Fits an empirical Bayes matrix factorization (see \strong{Details} for a
#'   description of the model). The resulting fit is referred to as a "flash"
#'   object (short for Factors and Loadings using Adaptive Shrinkage). Two
#'   interfaces are provided. The \code{flash} function provides a simple
#'   interface that allows a flash object to be fit in a single pass, while
#'   \code{flash.xxx} functions are pipeable functions that allow more
#'   complex flash objects to be fit incrementally. See the vignettes and the
#'   examples below for usage. See \strong{Details} for a discussion of the
#'   terminology used throughout the documentation.
#'
#' If \eqn{Y} is an \eqn{n \times p} data matrix, then the rank-one
#'   empirical Bayes matrix factorization model is:
#' \deqn{Y = \ell^{(1)} (\ell^{(2)})^T + E,} where \eqn{\ell^{(1)}} is an
#'   \eqn{n}-vector of \strong{row loadings}, \eqn{\ell^{(2)}} is a
#'   \eqn{p}-vector of \strong{column loadings}, and \eqn{E} is an
#'   \eqn{n \times p} matrix of \strong{residuals}. The following priors are
#'   assumed:
#' \deqn{E \sim N(0, s_{ij}^2)}
#' \deqn{\ell^{(1)} \sim g^{(1)}}
#' \deqn{\ell^{(2)} \sim g^{(2)}}
#' The residual variance parameters \eqn{s_{ij}^2} are constrained to have
#'   a simple structure (for example, the "constant" variance structure
#'   \eqn{s_{ij}^2 = s^2} for some \eqn{s^2} and for all \eqn{i}, \eqn{j})
#'   and fit via maximum likelihood.
#' The functions \eqn{g^{(1)}} and \eqn{g^{(2)}} are constrained to belong to
#'   some families of priors \eqn{G^{(1)}} and \eqn{G^{(1)}} (e.g., the family
#'   of point-normal distributions) and fit via variational approximation.
#'
#' The general rank-\eqn{K} empirical Bayes tensor factorization model is:
#'   \deqn{Y_{ijm} = \sum_k L^{(1)}_{ik} L^{(2)}_{jk} L^{(3)}_{mk} + E_{ijm},}
#'   where the \eqn{L^{(n)}}s are matrices of \strong{mode-n loadings}.
#' In the terminology used throughout this package, the \eqn{k}th column of
#'   \eqn{L^{(n)}} gives the mode-\eqn{n} loadings for the \eqn{k}th
#'   factor. That is, a \strong{factor} is a rank-one component of a matrix
#'   decomposition, formed by the outer product of vectors
#'   \eqn{L^{(1)}_{\cdot k}}, \eqn{L^{(2)}_{\cdot k}}, and
#'   \eqn{L^{(3)}_{\cdot k}}.
#'
#' Separate priors \eqn{g^{(n)}_k} are estimated for each mode and each factor,
#'   and different prior families \eqn{G^{(n)}_k} may be used for different
#'   modes and different factors. In general:
#' \deqn{E \sim N(0, s_{ijm}^2)}
#' \deqn{\L^{(n}_{\cdot k} \sim g^{(n)}_k \in G^{(n)}_k}
#'
#' @param data The observations. Can be a matrix, sparse matrix of class
#'   \code{Matrix}, or three-dimensional array. Can also be a low-rank matrix
#'   representation as returned by, for example, \code{\link{svd}},
#'   \code{\link[irlba]{irlba}}, \code{\link[rsvd]{rsvd}}, or
#'   \code{\link[softImpute]{softImpute}}.
#'
#' @param S The standard errors. Can be a matrix, scalar (if standard errors
#'   are the same for all observations), or vector (if, for example, \code{S}
#'   only varies across columns and is constant within any given row). If
#'   \code{NULL}, all residual variance will be estimated.
#'
#' @param prior.family Indicates the family of distributions that the priors on
#'   the loadings are
#'   assumed to belong to. Can be a list of length 1 or length \eqn{N}, where
#'   \eqn{N} is the number of modes (\eqn{N = 2} for matrices; \eqn{N = 3} for
#'   tensors). Each list element must be a prior family defined by one of the
#'   convenience functions \code{\link{prior.normal}},
#'   \code{\link{prior.point.normal}},
#'   \code{\link{prior.point.laplace}}, \code{\link{prior.nonzero.mode}},
#'   \code{\link{prior.scale.normal.mix}}, \code{\link{prior.unimodal}},
#'   \code{\link{prior.symmetric.unimodal}},
#'   \code{\link{prior.nonnegative}}, or \code{\link{prior.nonpositive}},
#'   or a custom prior type of a similar form (see \code{\link{prior.normal}}
#'   for details).
#'
#'   \code{prior.family} can also be a list of lists, in which case the first
#'   list specifies the family or families for the first factor, the second
#'   specifies the family or families for the second factor, and so on. The
#'   last list element is then re-used as often as necessary. See
#'   \code{\link{flash}} for examples.
#'
#' @param var.type Describes the structure of the estimated residual variance.
#'   Can be \code{NULL}, \code{0}, or a vector. If \code{NULL}, then
#'   \code{S} accounts for all residual variance. Otherwise, a rank-one
#'   variance structure will be estimated (and added to any variance specified
#'   by \code{S}). \code{var.type} then gives the modes along which the
#'   residual variance is permitted to vary. For example, \code{var.type = 1}
#'   estimates row-specific residual variances, while \code{var.type = c(1, 2)}
#'   optimizes over all rank-one matrices. If \code{var.type = 0}, then the
#'   residual variance is assumed to be constant across all observations.
#'
#'   Note that is usually faster to set \code{S = NULL} and to let \code{flash}
#'   estimate the residual variances. Further, it is much faster to only allow
#'   the residual variances to vary along a single mode. For example, estimating
#'   row-specific residual variances (\code{var.type = 1}) is much easier than
#'   estimating an arbitrary rank-one matrix of variances
#'   (\code{var.type = c(1, 2)}).
#'
#' @param greedy.Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by \code{flash}, since
#'   factors are only added as long as they increase the variational lower
#'   bound on the log likelihood for the model.
#'
#' @param backfit A "greedy" fit is performed by adding as many as
#'   \code{greedy.Kmax} factors, optimizing each newly added factor in one go
#'   without returning to optimize previously added factors. When
#'   \code{backfit = TRUE}, \code{flash} will additionally perform a final
#'   "backfit" where all factors are cyclically updated until convergence.
#'
#' @param nullcheck If \code{nullcheck = TRUE}, then \code{flash} will check
#'   that each factor in the final flash object improves the overall fit. Any
#'   factor that fails the check will be removed.
#'
#' @param verbose.lvl When and how to display progress updates. Set to
#'   \code{0} for none, \code{1} for updates after a factor is added or a
#'   backfit is completed, \code{2} for additional notifications about the
#'   variational lower bound, and \code{3} for updates after every iteration.
#'   Set to \code{-1} to output a single tab-delimited table of values.
#'
#' @return A \code{flash} object. Contains elements:
#'   \describe{
#'     \item{\code{n.factors}}{The total number of factors in the fitted
#'       model.}
#'     \item{\code{pve}}{The proportion of variance explained by each factor.}
#'     \item{\code{elbo}}{The variational lower bound achieved by the
#'       fitted model.}
#'     \item{\code{loadings.scale, loadings.pm}}{Posterior means for loadings.
#'       Since the model is not identifiable,
#'       each column of loadings is \eqn{L2}-normalized. The
#'       normalization constant is given by \code{loadings.scale}. Thus, for
#'       matrices, fitted values can be calculated as \code{f$loadings.pm[[1]]
#'       \%*\% diag(f$loadings.scale) \%*\% t(f$loadings.pm[[2]])} (or, more
#'       simply, as \code{fitted(f)}).}
#'     \item{\code{loadings.psd}}{Posterior standard deviations for loadings.}
#'     \item{\code{loadings.lfsr}}{Local false sign rates for loadings.}
#'     \item{\code{residuals.sd}}{Estimated residual standard deviations (these
#'       include any variance component given as an argument to \code{S}).}
#'     \item{\code{fitted.g}}{The fitted priors for each mode and factor.}
#'     \item{\code{sampler}}{A function that takes a single argument
#'       \code{nsamp} and returns \code{nsamp} samples from the posterior
#'       distribution of the (non-normalized) loadings.}
#'     \item{\code{flash.fit}}{A \code{flash.fit} object. Used by
#'       \code{flash} when fitting is not performed all at once, but
#'       incrementally via calls to various \code{flash.xxx} functions.}
#'   }
#'
#' @examples
#' data(gtex)
#'
#' # Fit up to 10 factors and backfit.
#' fl <- flash(gtex, greedy.Kmax = 10L, backfit = TRUE)
#'
#' # This is equivalent to the series of calls:
#' fl <- flash.init(gtex) %>%
#'   flash.add.greedy(Kmax = 10L) %>%
#'   flash.backfit() %>%
#'   flash.nullcheck()
#'
#' # Fit a unimodal distribution with mode zero and nonnegative support to
#' #   each set of row loadings and a scale mixture of normals with mean zero
#' #   to each set of column loadings.
#' fl <- flash(gtex,
#'             prior.family = c(prior.nonnegative(), prior.normal.scale.mix()),
#'             greedy.Kmax = 5)
#'
#' # Fit point-normal priors with nonzero means for the first factor and scale
#' #   mixtures of normals for every subsequent factor.
#' fl <- flash(gtex,
#'             prior.family = list(prior.nonzero.mode(), prior.normal.scale.mix()),
#'             greedy.Kmax = 5)
#'
#' # Fit a "Kronecker" (rank-one) variance structure (this can be slow).
#' fl <- flash(gtex, var.type = c(1, 2), greedy.Kmax = 5L)
#'
#' @seealso \code{\link{flash.init}}, \code{\link{flash.add.greedy}},
#'   \code{\link{flash.backfit}}, and \code{\link{flash.nullcheck}}. For more
#'   advanced functionality, see \code{\link{flash.init.factors}},
#'   \code{\link{flash.fix.loadings}}, \code{\link{flash.set.factors.to.zero}},
#'   \code{\link{flash.remove.factors}}, and \code{\link{flash.set.verbose}}.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
flash <- function(data,
                  S = NULL,
                  prior.family = prior.point.normal(),
                  var.type = 0L,
                  greedy.Kmax = 50L,
                  backfit = FALSE,
                  nullcheck = TRUE,
                  verbose.lvl = 1L) {
  fl <- flash.init(data, S = S, var.type = var.type)

  fl <- flash.add.greedy(fl, Kmax = greedy.Kmax, prior.family = prior.family,
                         verbose.lvl = verbose.lvl)

  if (backfit) {
    fl <- flash.backfit(fl, verbose.lvl = verbose.lvl)
  }

  if (nullcheck) {
    fl <- flash.nullcheck(fl, remove = TRUE, verbose.lvl = verbose.lvl)
  }

  return(fl)
}
