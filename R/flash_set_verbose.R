#' TODO: document
#'
#' @export
#'
flash.set.verbose <- function(flash,
                              verbose,
                              fns = NULL,
                              colnames = NULL,
                              colwidths = NULL) {
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

    fns       <- look.up.verbose.fns(verbose, data.dim)
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

  flash <- set.verbose.options(flash, lvl, fns, colnames, colwidths)

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
