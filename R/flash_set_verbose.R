#' Set verbose output
#'
#' Sets the default verbosity level and verbose output for a flash object.
#'
#' @param flash A \code{flash} or \code{flash.fit} object.
#'
#' @param verbose When and how to display progress updates. Set to
#'   \code{0} for none, \code{1} for updates after a factor is added or a
#'   backfit is completed, \code{2} for additional notifications about the
#'   variational lower bound, and \code{3} for updates after every iteration.
#'   Set to \code{-1} to output a single tab-delimited table of values.
#'
#' When \code{verbose = 3}, the improvement in the ELBO and the largest
#'   change in the loadings is shown after each iteration. If \code{verbose}
#'   is a character string rather than an integer, then a different set of
#'   columns will be shown. The character string should be a set of letters
#'   separated by spaces. Each letter may be followed by a number to specify
#'   a mode (e.g., rows or columns). Current options are: the (O)bjective
#'   after each iteration; the (D)ifference in the objective after each
#'   iteration; the maximum change in the mode-n (L)oadings after each
#'   iteration; the index corresponding to (W)hich mode-n loadings experienced
#'   the largest change; and the (S)parsity of the mode-n loadings after each
#'   iteration. See below for examples.
#'
#' Note that output is completely customizable via parameters \code{fns},
#'   \code{colnames}, and \code{colwidths}. If these are used, then
#'   \code{verbose} must take an integer argument (and the specified columns
#'   will only display when \code{verbose} is either -1 or 3).
#'
#' @param fns A vector of functions. Used to calculate values to output after
#'   each factor update when \code{verbose} is either -1 or 3.
#'
#' @param colnames A vector of column names, one for each function in
#'   \code{fns}.
#'
#' @param colwidths A vector of column widths.
#'
#' @examples
#' data(gtex)
#'
#' # Suppress all verbose output.
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose(0) %>%
#'   flash.add.greedy(Kmax = 5L)
#'
#' # Set verbose output using a character string.
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose("D L2 W2 S1") %>%
#'   flash.add.greedy(Kmax = 5L)
#'
#' # Set custom verbose output.
#' fns <- c(flashier:::calc.obj.diff,
#'          function(new, old, k) {
#'            flashier:::calc.max.abs.chg.EF(new, old, k, n = 2)
#'          })
#' colnames <- c("ELBO Diff", "Max Chg (Tiss)")
#' colwidths <- c(12, 18)
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose(3L, fns = fns, colnames = colnames, colwidths = colwidths) %>%
#'   flash.add.greedy(Kmax = 3L)
#'
#' # Default output can be set and then turned off as needed.
#' fl <- flash.init(gtex) %>%
#'   flash.set.verbose("D L2 W2 S1") %>%
#'   flash.add.greedy(Kmax = 5L) %>%
#'   flash.backfit(verbose.lvl = 1) %>%
#'   flash.add.greedy(Kmax = 1L)
#'
#' @export
#'
flash.set.verbose <- function(flash,
                              verbose,
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
    lvl <- verbose

    if (length(fns) == 0 && verbose == -1) {
      # Default output columns for verbose.lvl = -1.
      fns       <- c(get.new.obj, calc.obj.diff, calc.max.chg.EF)
      colnames  <- c("Obj", "Obj.diff", "Max.chg")
      colwidths <- c(14, 12, 12)
    } else if (length(fns) == 0) {
      # Default output columns for verbose.lvl = 3.
      fns       <- c(calc.obj.diff, calc.max.chg.EF)
      colnames  <- c("Obj Diff", "Max Chg")
      colwidths <- c(12, 12)
    }
  }

  fit <- set.verbose.options(fit, lvl, fns, colnames, colwidths)
  flash <- set.fit(flash, fit)

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
        return(calc.obj.diff)
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
