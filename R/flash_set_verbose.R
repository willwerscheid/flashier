#' Set verbose output
#'
#' Used in a \code{\link{flash}} pipeline to set the output that will be printed
#'   after each greedy or backfitting iteration.
#'
#' @details Function \code{flash_set_verbose} can be used to customize
#'   the output that is printed to console while fitting a \code{flash} object.
#'   After each greedy or backfitting iteration (see, respectively,
#'   \code{\link{flash_greedy}} and \code{\link{flash_backfit}}), each
#'   function in argument \code{fns} is successively evaluated and the
#'   result is printed to console in a table with column names defined by
#'   argument \code{colnames} and column widths defined by argument
#'   \code{colwidths}.
#'
#'   Each function in \code{fns} must accept exactly three parameters,
#'   \code{curr}, \code{prev}, and \code{k}: \code{curr} refers to the
#'   \code{\link{flash_fit}} object from the current iteration; \code{prev},
#'   to the \code{flash_fit} object from the previous iteration;
#'   and, if the iteration is a sequential backfitting iteration (that is, a
#'   \code{\link{flash_backfit}} iteration with argument
#'   \code{extrapolate = FALSE}), \code{k} identifies the factor/loadings pair
#'   that is currently being updated (in all other cases, \code{k} is
#'   \code{NULL}). Package \code{flashier} provides a number of functions that
#'   may be used to customize output: see
#'   \code{\link{flash_verbose_elbo}},
#'   \code{\link{flash_verbose_elbo_diff}},
#'   \code{\link{flash_verbose_max_chg}},
#'   \code{\link{flash_verbose_max_chg_L}}, and
#'   \code{\link{flash_verbose_max_chg_F}}. Custom functions may also be
#'   defined. They might inspect the current \code{flash_fit} object
#'   via argument \code{curr}; compare the fit in \code{curr} to the fit from the
#'   previous iteration (provided by argument \code{prev}); or
#'   ignore both \code{flash_fit} objects entirely (for example, to
#'   track progress over time, one might simply call \code{\link{Sys.time}}).
#'   To facilitate working with \code{flash_fit} objects, package
#'   \code{flashier} provides a number of accessors, which are enumerated in
#'   the documentation for object \code{\link{flash_fit}}. Custom functions
#'   should return a character string that contains the output exactly as it is
#'   to displayed; see \strong{Examples} below.
#'
#' @param flash A \code{flash} or \code{flash_fit} object.
#'
#' @param verbose When and how to display progress updates. Set to \code{0}
#'   for no updates; \code{1} for updates after a "greedy" factor is added or
#'   a backfit is completed; \code{2} for additional notifications about the
#'   variational lower bound (ELBO); and \code{3} for updates after every
#'   iteration. By default, per-iteration update information includes the
#'   change in ELBO and the maximum (absolute) change over all L2-normalized
#'   loadings \eqn{\ell_1, \ldots, \ell_K} and factors \eqn{f_1, \ldots, f_K}.
#'   Update information is customizable via parameters \code{fns},
#'   \code{colnames}, and \code{colwidths}.
#'
#' A single tab-delimited table of values may also be output using
#'   option \code{verbose = -1}. This format is especially convenient for
#'   downstream analysis of the fitting history. For example, it may be used
#'   to plot the value of the ELBO after each iteration (see the "Advanced
#'   Flashier" vignette for an illustration).
#'
#' @param fns A vector of functions. Used to calculate values to display
#'   after each greedy/backfit iteration when \code{verbose} is either -1 or 3
#'   (see Details below). Ignored for other values of \code{verbose} (0, 1, or 2).
#'
#' @param colnames A vector of column names, one for each function in
#'   \code{fns}.
#'
#' @param colwidths A vector of column widths, one for each function in
#'   \code{fns}.
#'
#' @return The \code{\link{flash}} object from argument \code{flash}, with the
#'   new verbose settings reflected in updates to the "internal"
#'   \code{flash_fit} object. These settings will persist across
#'   all subsequent calls to \code{flash_xxx} functions until they are modified
#'   by another call to \code{flash_set_verbose}.
#'
#' @examples
#' # Suppress all verbose output.
#' fl <- flash_init(gtex) |>
#'   flash_set_verbose(0) |>
#'   flash_greedy(Kmax = 5)
#'
#' # Set custom verbose output.
#' sparsity_F <- function(curr, prev, k) {
#'   g_F <- flash_fit_get_g(curr, n = 2)
#'   g_F_pi0 <- g_F$pi[1] # Mixture weight of the "null" component.
#'   return(g_F_pi0)
#' }
#' verbose_fns <- c(flash_verbose_elbo, flash_verbose_max_chg_F, sparsity_F)
#' colnames <- c("ELBO", "Max Chg (Tiss)", "Sparsity (Tiss)")
#' colwidths <- c(12, 18, 18)
#' fl <- flash_init(gtex) |>
#'   flash_set_verbose(
#'     verbose = 3,
#'     fns = verbose_fns,
#'     colnames = colnames,
#'     colwidths = colwidths
#'   ) |>
#'   flash_greedy(Kmax = 3)
#'
#' # Output can be changed as needed.
#' fl <- flash_init(gtex) |>
#'   flash_set_verbose(verbose = 1) |>
#'   flash_greedy(Kmax = 5L) |>
#'   flash_backfit(verbose = 3) |>
#'   flash_greedy(Kmax = 1L)
#'
#' @export
#'
flash_set_verbose <- function(flash,
                              verbose = 1L,
                              fns = NULL,
                              colnames = NULL,
                              colwidths = NULL) {
  fit <- get.fit(flash)

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
      fns       <- c(flash_verbose_elbo, flash_verbose_elbo_diff, flash_verbose_max_chg)
      colnames  <- c("ELBO", "ELBO.diff", "LF.max.chg")
      colwidths <- c(14, 12, 12)
    } else if (length(fns) == 0) {
      # Default output columns for verbose.lvl = 3.
      fns       <- c(flash_verbose_elbo_diff, flash_verbose_max_chg)
      colnames  <- c("ELBO Diff", "LF Max Chg")
      colwidths <- c(12, 12)
    } else if (missing(verbose)) {
      verbose <- 3
    }

    lvl <- verbose
  }

  fit <- set.verbose.options(fit, lvl, fns, colnames, colwidths)
  flash <- set.fit(flash, fit)

  if (lvl == -1) {
    print_tab.delim.table.header(colnames)
  }

  return(flash)
}

look.up.verbose.fns <- function(verbose, data.dim) {
  fns <- lapply(verbose, function(symbol) {
    # There is the option to pass in a character string as shorthand; this
    #   option has not yet been documented except in the advanced vignette.
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
        return(flash_verbose_elbo)
      } else { # if chars[[1]] == "D"
        return(flash_verbose_elbo_diff)
      }
    } else if (chars[[1]] == "L") {
      return(function(curr, prev, k) calc.max.chg.EF(curr, prev, k, n))
    } else if (chars[[1]] == "W") {
      return(function(curr, prev, k) which.max.chg.EF(curr, prev, k, n))
    } else if (chars[[1]] == "S") {
      if (is.null(n))
        stop("Dimension must be specified for verbose sparsity output.")
      return(function(curr, prev, k) get.sparsity(curr, prev, k, n))
    } else if (chars[[1]] == "E") {
      return(function(curr, prev, k) get.exclusion.count(curr, prev, k, n))
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
         "flash_set_verbose with verbose = -1. See ?flash_set_verbose ",
         "for usage.")
  } else {
    verbose.lvl <- verbose
  }

  return(verbose.lvl)
}
