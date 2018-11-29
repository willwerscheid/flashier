# Level 1 announcements (main loop) -------------------------------------------

announce.flash.init <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    message("Initializing flash object...")
}

announce.add.factor <- function(verbose.lvl, k) {
  if (verbose.lvl > 0)
    message("Adding factor ", k, " to flash object...")
}

report.greedy.obj.decrease <- function(verbose.lvl, obj.diff) {
  if (verbose.lvl > 0)
    message("An iteration decreased the objective by ",
            formatC(-obj.diff, format = "e", digits = 2),
            ". Try backfitting with warmstarts.")
}

report.add.factor.result <- function(verbose.lvl, greedy.complete, obj) {
  if (verbose.lvl > 0 && greedy.complete) {
    message("Factor doesn't increase objective and won't be added.")
  }
  if (verbose.lvl > 1 && !greedy.complete) {
    message("  Factor successfully added. Objective: ",
            formatC(obj, format = "f", digits = 3))
  }
}

announce.backfit <- function(verbose.lvl, n.factors) {
  if (verbose.lvl > 0)
    message("Backfitting ", n.factors, " factors...")
}

report.backfit.obj.decrease <- function(verbose.lvl, obj.diff, k) {
  if (verbose.lvl > 0)
    message("An update to factor ", k, " decreased the objective by ",
            formatC(-obj.diff, format = "e", digits = 3),
            ". Try using warmstarts?")
}

announce.nullchk <- function(verbose.lvl, n.factors) {
  if (verbose.lvl > 0 && n.factors > 0)
    message("Nullchecking ", n.factors, " factors...")
}

report.nullchk.failure <- function(verbose.lvl, obj.diff, k) {
  if (verbose.lvl > 0) {
    if (obj.diff > 0) {
      message("Factor ", k, " removed, increasing objective by ",
              formatC(obj.diff, format = "e", digits = 3), ".")
    } else if (obj.diff == 0) {
      message("Factor ", k, " removed with no change to objective.")
    }
  }
}

announce.wrapup <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    message("Wrapping up...")
}

report.completion <- function(verbose.lvl) {
  if (verbose.lvl > 0)
    message("Done.")
}

# Level 2 announcements -------------------------------------------------------

announce.factor.init <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    message("  Initializing factor...")
}

announce.factor.opt <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    message("  Optimizing factor...")
}

report.backfit.complete <- function(verbose.lvl, obj) {
  if (verbose.lvl > 1)
    message("  Backfit complete. Objective: ",
            formatC(obj, format = "f", digits = 3))
}

report.nullchk.success <- function(verbose.lvl) {
  if (verbose.lvl > 1)
    message("  No factor can be removed without decreasing the objective.")
}

# Optimization details (level 3) ----------------------------------------------

report.tol.setting <- function(verbose.lvl, tol) {
  if (verbose.lvl > 2)
    message("Convergence tolerance set to ",
            formatC(tol, format = "e", digits = 2), ".")
}

print.table.header <- function(verbose.lvl, colnames, colwidths, backfit) {
  if (verbose.lvl > 2) {
    header.string <- sprintf("%13s", "Iteration")
    if (backfit)
      header.string <- paste0(header.string, sprintf("%8s", "Factor"))
    for (col in 1:length(colnames)) {
      width.string <- paste0("%", as.character(colwidths[col]), "s")
      header.string <- paste0(header.string,
                              sprintf(width.string, colnames[col]))
    }
    message(header.string)
  }
}

print.tab.delim.table.header <- function(colnames) {
  header.string <- "Type\tFactor\tIter"
  for (name in colnames) {
    header.string <- paste0(header.string, "\t", name)
  }
  header.string <- paste0(header.string, "\n")
  cat(header.string)
}

print.table.entry <- function(verbose.lvl, colwidths, iter, info, k, backfit) {
  if (verbose.lvl > 2) {
    table.entry <- sprintf("%13d", iter)
    if (backfit)
      table.entry <- paste0(table.entry, sprintf("%8d", k))
    for (col in 1:length(colwidths)) {
      width.string <- paste0("%", as.character(colwidths[col]), "s")
      if (is.finite(info[col]) && round(info[col]) == info[col]) {
        format.info <- formatC(info[col], format = "d")
      } else {
        format.info <- formatC(info[col], format = "e", digits = 2)
      }
      table.entry <- paste0(table.entry, sprintf(width.string, format.info))
    }
    message(table.entry)
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
