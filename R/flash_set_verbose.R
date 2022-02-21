#' Set verbose output
#'
#' Determines the output that will be displayed when fitting a \code{flash}
#'   object.
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param verbose When and how to display progress updates. Set to \code{0}
#'   for no updates; \code{1} for updates after a "greedy" factor is added or
#'   a backfit is completed; \code{2} for additional notifications about the
#'   variational lower bound (ELBO); and \code{3} for updates after every
#'   iteration. By default, per-iteration update information includes the
#'   change in ELBO and the maximum (absolute) change over all L2-normalized
#'   loadings \eqn{\ell_1, \ldots, \ell_K} and factors \eqn{f_1, \ldots, f_K}.
#'   Update information is customizable via parameters \code{disp.fns},
#'   \code{colnames}, and \code{colwidths}.
#'
#' A single tab-delimited table of values may also be output using
#'   option \code{verbose = -1}. This format is especially convenient for
#'   downstream analysis of the fitting history (for example, it may be used
#'   to plot the value of the ELBO after each iteration).
#'
#' @param disp.fns A vector of functions. Used to calculate values to display
#'   after each greedy/backfit iteration when \code{verbose} is either -1 or 3.
#'   Options include \code{\link{display.elbo}},
#'   \code{\link{display.elbo.diff}}, \code{\link{display.max.chg}},
#'   \code{\link{display.L.max.chg}}, and \code{\link{display.F.max.chg}}.
#'   Custom functions may also be used. They should accept three parameters,
#'   \code{new}, \code{old}, and \code{k}, where \code{new} refers to the
#'   \code{\link{flash.fit}} object from the current iteration, \code{old}
#'   refers to the \code{flash.fit} object from the previous iteration,
#'   and \code{k} identifies the factor/loadings pair that is currently
#'   being updated during sequential backfits (that is, in calls to function
#'   \code{\link{flash.backfit}} where \code{extrapolate = FALSE}). See below
#'   for an example.
#'
#' @param colnames A vector of column names, one for each function in
#'   \code{disp.fns}.
#'
#' @param colwidths A vector of column widths, one for each function in
#'   \code{disp.fns}.
#'
#' @return A \code{\link{flash}} object.
#'
#' @seealso \code{\link{display.elbo}}, \code{\link{display.elbo.diff}},
#'   \code{\link{display.max.chg}}, \code{\link{display.L.max.chg}},
#'   \code{\link{display.F.max.chg}}
#'
#' @examples
#' data(gtex)
#'
#' # Suppress all verbose output.
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose(0) %>%
#'   flash.add.greedy(Kmax = 5L)
#'
#' # Set custom verbose output.
#' sparsity.F <- function(new, old, k) {
#'   g.F <- ff.g(new, k, n = 2)
#'   g.F.pi0 <- g.F$pi[1] # Mixture weight of the "null" component.
#'   return(g.F.pi0)
#' }
#' disp.fns <- c(display.elbo, display.F.max.chg, sparsity.F)
#' colnames <- c("ELBO", "Max Chg (Tiss)", "Sparsity (Tiss)")
#' colwidths <- c(12, 18, 18)
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose(
#'     3L,
#'     disp.fns = disp.fns,
#'     colnames = colnames,
#'     colwidths = colwidths
#'   ) %>%
#'   flash.add.greedy(Kmax = 3L)
#'
#' # Output can be changed as needed.
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose(verbose = 1) %>%
#'   flash.add.greedy(Kmax = 5L) %>%
#'   flash.backfit(verbose = 3) %>%
#'   flash.add.greedy(Kmax = 1L)
#'
#' @export
#'
flash.set.verbose <- function(flash,
                              verbose = 1L,
                              disp.fns = NULL,
                              colnames = NULL,
                              colwidths = NULL) {
  fit <- get.fit(flash)
  fns <- disp.fns

  if (!all(length(fns) == c(length(colnames), length(colwidths)))) {
    stop("Arguments to fns, colnames, and colwidths must all have the same ",
         "length.")
  }

  if (is.character(verbose)) {
    if(length(fns) > 0) {
      stop("If verbose output is specified via a character argument to ",
           "verbose, then fns, colnames, and colwidths must be null.")
    }

    verbose <- unlist(strsplit(toupper(verbose), "[ .,/]"))
    verbose <- verbose[verbose != ""]

    lvl <- 3

    fns       <- look.up.verbose.fns(verbose, get.dim(fit))
    colnames  <- look.up.verbose.colnames(verbose)
    colwidths <- look.up.verbose.colwidths(verbose)
  } else {
    if (length(fns) == 0 && verbose == -1) {
      # Default output columns for verbose.lvl = -1.
      fns       <- c(display.elbo, display.elbo.diff, display.max.chg)
      colnames  <- c("ELBO", "ELBO.diff", "Max.chg")
      colwidths <- c(14, 12, 12)
    } else if (length(fns) == 0) {
      # Default output columns for verbose.lvl = 3.
      fns       <- c(display.elbo.diff, display.max.chg)
      colnames  <- c("ELBO Diff", "Max Chg")
      colwidths <- c(12, 12)
    } else if (missing(verbose)) {
      verbose <- 3
    }

    lvl <- verbose
  }

  fit <- set.verbose.options(fit, lvl, fns, colnames, colwidths)
  flash <- set.fit(flash, fit)

  if (lvl == -1) {
    print.tab.delim.table.header(colnames)
  }

  return(flash)
}

look.up.verbose.fns <- function(verbose, data.dim) {
  fns <- lapply(verbose, function(symbol) {
    chars <- unlist(strsplit(symbol, ""))

    if (length(chars) > 2) {
      stop("Unable to parse verbose output string.")
    } else if (length(chars) == 2) {
      n <- as.integer(chars[[2]])
      must.be.integer(n, lower = 1, upper = data.dim)
    } else {
      n <- NULL
    }

    if (chars[[1]] %in% c("O", "D")) {
      if (!is.null(n))
        warning("Dimension ignored for verbose objective output.")
      if (chars[[1]] == "O") {
        return(display.obj)
      } else { # if chars[[1]] == "D"
        return(display.obj.diff)
      }
    } else if (chars[[1]] == "L") {
      return(function(new, old, k) calc.max.chg.EF(new, old, k, n))
    } else if (chars[[1]] == "W") {
      return(function(new, old, k) which.max.chg.EF(new, old, k, n))
    } else if (chars[[1]] == "S") {
      if (is.null(n))
        stop("Dimension must be specified for verbose sparsity output.")
      return(function(new, old, k) get.sparsity(new, old, k, n))
    } else if (chars[[1]] == "E") {
      return(function(new, old, k) get.exclusion.count(new, old, k, n))
    }

    stop("Unrecognized verbose output character.")
  })

  return(fns)
}

look.up.verbose.colnames <- function(verbose) {
  names <- lapply(verbose, function(symbol) {
    chars <- unlist(strsplit(symbol, ""))

    if (chars[[1]] == "O") {
      return("Objective")
    } else if (chars[[1]] == "D") {
      return("Obj Diff")
    } else if (chars[[1]] == "L") {
      name <- "Max Chg"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    } else if (chars[[1]] == "W") {
      name <- "Whch"
      if (length(chars) > 1)
        name <- paste(name, chars[[2]])
      return(name)
    } else if (chars[[1]] == "S") {
      return(paste("Sparsity", chars[[2]]))
    } else if (chars[[1]] == "E") {
      return(paste("Excl", chars[[2]]))
    }

    stop("Unrecognized verbose output character.")
  })

  return(unlist(names))
}

look.up.verbose.colwidths <- function(verbose) {
  widths <- lapply(verbose, function(symbol) {
    char <- substring(symbol, 1, 1)

    if (char == "O") {
      return(16)
    } else if (char %in% c("D", "L", "S")) {
      return(13)
    } else if (char %in% c("W", "E")) {
      return(9)
    }

    stop("Unrecognized verbose output character.")
  })

  return(unlist(widths))
}

handle.verbose.param <- function(verbose, flash) {
  if (is.null(verbose)) {
    verbose.lvl <- get.verbose.lvl(flash)
  } else if (verbose == -1) {
    stop("To output a tab-delimited table of values, use function ",
         "flash.set.verbose with verbose = -1. See ?flash.set.verbose ",
         "for usage.")
  } else {
    verbose.lvl <- verbose
  }

  return(verbose.lvl)
}
