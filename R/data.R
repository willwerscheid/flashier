#' @name gtex
#'
#' @title GTEx data
#'
#' @docType data
#'
#' @description Derived from data made available by the Genotype Tissue
#'   Expression (GTEx) project (Lonsdale et al. 2013),
#'   which provides \eqn{z}-scores for assessing the significance of effects of
#'   genetic variants (single nucleotide polymorphisms, or SNPs) on gene
#'   expression across 44 human tissues. To reduce the data to a more manageable
#'   size, Urbut et al. (2019) chose the "top" SNP for
#'   each gene --- that is, the SNP associated with the largest (absolute)
#'   \eqn{z}-score over all 44 tissues. This yields a \eqn{16,069 \times 44} matrix
#'   of \eqn{z}-scores, with rows corresponding to SNP-gene pairs and columns
#'   corresponding to tissues. The dataset included here
#'   is further subsampled down to 1000 rows.
#'
#' @format
#' \code{gtex} is a matrix with 1000 rows and 44 columns, with rows
#'   corresponding to SNP-gene pairs and columns corresponding to tissues.
#'
#' @source <https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds>
#'
#' @references
#' Lonsdale et al. (2013).
#'   "The Genotype-Tissue Expression (GTEx) project." \emph{Nature Genetics}
#'   45(6), 580--585.
#'
#' Urbut, Wang, Carbonetto, and Stephens (2019).
#'   "Flexible statistical methods for estimating and testing effects in
#'   genomic studies with multiple conditions." \emph{Nature Genetics}
#'   51(1), 187--195.
#'
#' @keywords data
#'
#' @examples
#' data(gtex)
#' summary(gtex)
#'
NULL

#' @name gtex_colors
#'
#' @title Colors for plotting GTEx data
#'
#' @docType data
#'
#' @description A custom palette used by Wang and Stephens (2021) to plot an
#'   empirical Bayes matrix factorization of data from the GTEx project
#'   (of which the \code{\link{gtex}} data in package \strong{flashier} is a
#'   subsample).
#'   The palette is designed to link similar tissues together visually. For
#'   example, brain tissues all have the same color (yellow); arterial tissues
#'   are shades of pink or red; etc.
#'
#' @format
#' \code{gtex_colors} is a named vector of length 44, with names corresponding
#'   to tissues (columns) in the \code{\link{gtex}} dataset and values
#'   giving hexadecimal color codes.
#'
#' @source <https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt>
#'
#' @keywords data
#'
#' @references
#' Wei Wang and Matthew Stephens (2021).
#'   "Empirical Bayes matrix factorization." \emph{Journal of Machine Learning
#'   Research} 22, 1--40.
#'
#' @examples
#' fl <- flash(gtex, greedy_Kmax = 4)
#' plot(fl, pm_colors = gtex_colors)
NULL
