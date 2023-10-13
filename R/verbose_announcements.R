# Level 1 announcements (main loop) -------------------------------------------

announce.flash.init <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    cat("Initializing flash object...\n")
}

announce.add.factor <- function(verbose.lvl, k) {
  if (verbose.lvl > 0)
    cat("Adding factor", k, "to flash object...\n")
}

report.add.factor.result <- function(verbose.lvl, greedy.complete, obj) {
  if (verbose.lvl > 0 && greedy.complete) {
    cat("Factor doesn't significantly increase objective and won't be",
        "added.\n")
  }
  if (verbose.lvl > 1 && !greedy.complete) {
    cat("  Factor successfully added. Objective:",
        formatC(obj, format = "f", digits = 3), "\n")
  }
}

announce.backfit <- function(verbose.lvl, n.factors, tol) {
  if (verbose.lvl > 0)
    cat(paste0("Backfitting ", n.factors, " factors (tolerance: ",
               formatC(tol, format = "e", digits = 2), ")...\n"))
}

announce.no.backfit <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    cat("No factors have been added. Skipping backfit.\n")
}

report.backfit.progress <- function(verbose.lvl, tol) {
  if (verbose.lvl > 0 && verbose.lvl < 3) {
    cat(paste0("  Difference between iterations is within ",
               formatC(tol, format = "e", digits = 1),
               "...\n"))
  }
}

report.maxiter.reached <- function(verbose.lvl) {
  if (verbose.lvl > 0) {
    cat("  --Maximum number of iterations reached!\n")
  } else {
    warning("Maximum number of iterations reached.")
  }
}

report.timeout.reached <- function(verbose.lvl, t.diff) {
  if (verbose.lvl > 0) {
    cat("  --Fit timed out after", format(t.diff), "!\n")
  } else {
    warning("Fit timed out after ", format(t.diff), ".")
  }
}

report.timeout.no.greedy <- function(verbose.lvl) {
  if (verbose.lvl > 0) {
    cat("Fit timed out. Skipping greedy fit.\n")
  } else {
    warning("Fit timed out. Greedy fit skipped.")
  }
}

report.timeout.no.backfit <- function(verbose.lvl) {
  if (verbose.lvl > 0) {
    cat("Fit timed out. Skipping backfit.\n")
  } else {
    warning("Fit timed out. Backfit skipped.")
  }
}

report.backfit.obj.decrease <- function(verbose.lvl, obj.diff, k) {
  if (verbose.lvl > 0)
    cat("An update to factor ", k, " decreased the objective by ",
        formatC(-obj.diff, format = "e", digits = 3), ".\n", sep = "")
}

announce.nullchk <- function(verbose.lvl, n.factors, n.zero.factors) {
  if (verbose.lvl > 0 && n.factors > 0) {
    cat("Nullchecking", n.factors, "factors...\n")
    if (n.zero.factors == 1) {
      cat("  One factor is identically zero.\n")
    }
    if (n.zero.factors > 1) {
      cat(" ", n.zero.factors, "factors are identically zero.\n")
    }
  }
}

announce.no.nullchk <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    cat("No factors have been added. Skipping nullcheck.\n")
}

report.nullchk.failure <- function(verbose.lvl, obj.diff, k) {
  if (verbose.lvl > 0) {
    if (obj.diff > 0) {
      cat("Factor", k, "set to zero, increasing objective by ",
          formatC(obj.diff, format = "e", digits = 3), ".\n", sep = "")
    } else if (obj.diff == 0) {
      cat("Factor", k, "set to zero with no change to objective.\n")
    } else {
      cat("Factor", k, "set to zero, decreasing objective by ",
          formatC(-obj.diff, format = "e", digits = 3), ".\n", sep = "")
    }
  }
}

report.zero.factor <- function(verbose.lvl, k) {
  if (!is.null(k) && verbose.lvl > 0) {
    cat("  --Estimate of factor", k, "is numerically zero!\n")
  }
}

report.factor.removal <- function(verbose.lvl, n.factors) {
  if (verbose.lvl > 0) {
    if (n.factors == 1) {
      cat("  Removed one factor.\n")
    } else {
      cat("  Removed", n.factors, "factors.\n")
    }
  }
}

announce.wrapup <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    cat("Wrapping up...\n")
}

report.completion <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    cat("Done.\n")
}

# Level 2 announcements -------------------------------------------------------

announce.factor.init <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    cat("  Initializing factor...\n")
}

announce.factor.opt <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    cat("  Optimizing factor...\n")
}

report.backfit.complete <- function(verbose.lvl, obj) {
  if (verbose.lvl > 1)
    cat("  Backfit complete. Objective:",
        formatC(obj, format = "f", digits = 3), "\n")
}

report.nullchk.success <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    cat("  No factor can be removed without significantly decreasing the",
        "objective.\n")
}

# Optimization details (level 3) ----------------------------------------------

report.tol.setting <- function(verbose.lvl, tol) {
  if (verbose.lvl > 2)
    cat("Convergence tolerance set to ",
        formatC(tol, format = "e", digits = 2), ".\n", sep = "")
}

report.greedy.obj.decrease <- function(verbose.lvl, obj.diff) {
  if (verbose.lvl > 2)
    cat("An iteration decreased the objective by ",
        formatC(-obj.diff, format = "e", digits = 2), ".\n", sep="")
}

print_table.header <- function(verbose.lvl, colnames, colwidths, backfit) {
  if (verbose.lvl > 2) {
    header.string <- sprintf("%13s", "Iteration")
    if (backfit)
      header.string <- paste0(header.string, sprintf("%8s", "Factor"))
    for (col in 1:length(colnames)) {
      width.string <- paste0("%", as.character(colwidths[col]), "s")
      header.string <- paste0(header.string,
                              sprintf(width.string, colnames[col]))
    }
    header.string <- paste0(header.string, "\n")
    cat(header.string)
  }
}

print_tab.delim.table.header <- function(colnames) {
  header.string <- "Type\tFactor\tIter"
  for (name in colnames) {
    header.string <- paste0(header.string, "\t", name)
  }
  header.string <- paste0(header.string, "\n")
  cat(header.string)
}

print_table.entry <- function(verbose.lvl, colwidths, iter, info, k, backfit) {
  if (verbose.lvl > 2) {
    table.entry <- sprintf("%13d", iter)
    if (backfit)
      table.entry <- paste0(table.entry, sprintf("%8s", as.character(k)))
    for (col in 1:length(colwidths)) {
      width.string <- paste0("%", as.character(colwidths[col]), "s")
      if (is.numeric(info[[col]])
          && is.finite(info[[col]])
          && round(info[[col]]) == info[[col]]) {
        format.info <- formatC(info[[col]], format = "d")
      } else if (is.numeric(info[[col]])) {
        format.info <- formatC(info[[col]], format = "e", digits = 2)
      } else {
        format.info <- toString(info[[col]])
      }
      table.entry <- paste0(table.entry, sprintf(width.string, format.info))
    }
    cat(table.entry, "\n")
  } else if (verbose.lvl == -1) {
    if (backfit) {
      table.entry <- "backfit\t"
    } else {
      table.entry <- "greedy\t"
    }
    table.entry <- paste0(table.entry, k, "\t", iter)
    for (col in 1:length(colwidths)) {
      table.entry <- paste0(table.entry, "\t", info[col])
    }
    table.entry <- paste0(table.entry, "\n")
    cat(table.entry)
  }
}
