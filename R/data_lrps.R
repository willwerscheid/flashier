# Low-rank plus sparse matrix representations: X = UV' + S.
#
# Internally, an object of class "lrps" is a list of length two, with
# field "lowrank" containing an object of class "lowrank" (UV') and field
# "sparse" containing the matrix S, ideally but not necessarily of class
# "Matrix" (from the Matrix package).
#
#' @exportS3Method as.fullrank lrps
as.fullrank.lrps <- function(x, ...) {
  return(lowrank.expand(x[["lowrank"]]) + x[["sparse"]])
}

#' @exportS3Method get.data.dims lrps
get.data.dims.lrps <- function(x, ...) {
  return(dim(x[["sparse"]]))
}

#' @exportS3Method get.data.dimnames lrps
get.data.dimnames.lrps <- function(x, ...) {
  return(lapply(x[["lowrank"]], rownames))
}

#' @exportS3Method nmode.prod.vec lrps
nmode.prod.vec.lrps <- function(x, v, n, ...) {
  return(
    nmode.prod.vec(x[["lowrank"]], v, n) +
      nmode.prod.vec(x[["sparse"]], v, n)
  )
}

#' @exportS3Method premult.nmode.prod.r1 lrps
premult.nmode.prod.r1.lrps <- function(x, lr, r1, n, ...) {
  return(premult.nmode.prod.r1(x[["lowrank"]], lr, r1, n) +
           premult.nmode.prod.r1(x[["sparse"]], lr, r1, n))
}

#' @exportS3Method sq.nmode.prod.r1 lrps
sq.nmode.prod.r1.lrps <- function(x, r1, n, ...) {
  return(sq.nmode.prod.r1(x[["lowrank"]], r1, n) +
           sq.nmode.prod.r1(x[["sparse"]], r1, n) +
           2 * premult.nmode.prod.r1(x[["sparse"]], x[["lowrank"]], r1, n))
}

# Not currently used -----

# Only needed for parallel backfits:
#
nmode.prod.lowrank.lrps <- function(x, Y, n, ...) {
  return(nmode.prod.lowrank(x[["lowrank"]], Y, n) +
           nmode.prod.lowrank(x[["sparse"]], Y, n))
}
