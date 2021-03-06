% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flash_fix_loadings.R
\name{flash.fix.loadings}
\alias{flash.fix.loadings}
\title{Fix flash loadings}
\usage{
flash.fix.loadings(flash, kset, mode, is.fixed = TRUE)
}
\arguments{
\item{flash}{A \code{flash} or \code{flash.fit} object}

\item{kset}{A vector of integers indexing the factors whose loadings are to
be fixed (or unfixed).}

\item{mode}{The mode along which loadings are to be fixed (for matrix data,
\code{mode = 1} fixes row loadings and \code{mode = 2} fixes column
loadings). For any given factor, loadings may only be fixed along a single
mode.}

\item{is.fixed}{If \code{is.fixed = TRUE}, then all loadings along the
specified mode will be fixed. If only a subset of loadings are to be
fixed, then \code{is.fixed} should be an appropriately-sized vector or
matrix of values that can be coerced to logical. For example, if row
loadings for two factors are to be fixed, then \code{is.fixed} can be
a length-n vector or an n by 2 matrix (where n is the number of rows in
the data matrix). Finally, loadings can be "unfixed" by setting
\code{is.fixed = FALSE}.}
}
\description{
Fixes some or all loadings within one or more flash factors.
}
