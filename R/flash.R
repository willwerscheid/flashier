#' Empirical Bayes matrix/tensor factorization
#'
#' TODO: Describe model. Add reference.
#'
#' @param data The observations. Can be a matrix, sparse matrix of class
#'   \code{Matrix}, three-dimensional array, or \code{flash.data} object
#'   obtained from \code{\link{set.flash.data}}. Can also be a low-rank matrix
#'   representation as returned by, for example, \code{\link{svd}},
#'   \code{\link[irlba]{irlba}}, \code{\link[rsvd]{rsvd}}, or
#'   \code{\link[softImpute]{softImpute}}. Can be \code{NULL} if
#'   \code{init} is used.
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
#'   For example, the default \code{prior.family = prior.point.normal()} fits a
#'   (different) point-normal prior for each factor and each mode, while
#'   \code{prior.family = c(prior.nonnegative(), prior.scale.normal.mix())}
#'   fits a unimodal distribution with mode zero and nonnegative support to
#'   each set of row loadings and a scale mixture of normals with mean zero to
#'   each set of column loadings.
#'
#'   \code{prior.family} can also be a list of lists, in which case the first
#'   list specifies the family or families for the first factor, the second
#'   specifies the family or families for the second factor, and so on. The
#'   last list element is then re-used as often as necessary.
#'   For example, \code{prior.family = list(prior.nonzero.mode(),
#'   prior.scale.normal.mix())} will fit point-normal priors with nonzero
#'   means for the first factor and scale mixtures of normals for every
#'   subsequent factor.
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
#' @param backfit A "greedy" fit is performed by adding as many as
#'   \code{greedy.Kmax} factors, optimizing each newly added factor in one go
#'   without returning to optimize previously added factors. When
#'   \code{backfit = TRUE}, \code{flash} will additionally perform a final
#'   "backfit" where all factors are cyclically updated until convergence.
#'
#' @param greedy.Kmax The maximum number of factors to be added. This will not
#'   necessarily be the total number of factors added by \code{flash}, since
#'   factors are only added as long as they increase the variational lower
#'   bound on the log likelihood for the model. Fixed factors are not counted
#'   towards this limit.
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
#'     \item{\code{convergence.status}}{A character string indicating whether
#'       the fitting algorithm has converged.}
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
#'       incrementally via repeated calls to \code{flash} (with the
#'       intermediate \code{flash} or \code{flash.fit} objects given as
#'       arguments to \code{init}).}
#'   }
#'
#' @examples
#' # TODO
#'
#' @export
#'
flash <- function(data = NULL,
                  S = NULL,
                  prior.family = prior.point.normal(),
                  var.type = 0L,
                  backfit = FALSE,
                  greedy.Kmax = 50L,
                  verbose.lvl = 1L) {
  fl <- flash.init(data, var.type, S)

  fl <- flash.add.greedy(fl, greedy.Kmax, prior.family,
                         verbose.lvl = verbose.lvl,
                         output.lvl = 3L * !backfit)

  if (backfit) {
    fl <- flash.backfit(fl, verbose.lvl = verbose.lvl)
  }

  return(fl)
}
