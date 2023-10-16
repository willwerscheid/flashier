#' Plot method for flash objects
#'
#' Given a \code{\link{flash}} object, produces up to two figures: one showing
#'   the proportion of variance explained per factor/loadings pair, and one that
#'   plots posterior means for either factors or loadings (depending on the
#'   argument to parameter \code{pm_which}).
#'
#' @param x An object inheriting from class \code{flash}.
#'
#' @param include_scree Whether to include a figure ("scree plot") showing the
#'   proportion of variance explained by each factor/loadings pair.
#'
#' @param include_pm Whether to include a figure showing the posterior means for
#'   either loadings \eqn{L} or factors \eqn{F} (depending on the argument to
#'   parameter \code{pm_which}). One plot panel is produced for each
#'   factor/loadings pair \eqn{k}. If argument \code{pm_groups}
#'   is left unspecified, then bar plots will be produced, with each bar
#'   corresponding to a single value \eqn{\ell_{ik}} or \eqn{f_{jk}}.
#'   Otherwise, overlapping histograms will be
#'   produced, with each histogram corresponding to one of the groups
#'   specified by \code{pm_groups}.
#'
#' @param order_by_pve If \code{TRUE}, then the factor/loadings pairs will be
#'   re-ordered according to proportion of variance explained (from
#'   highest to lowest).
#'
#' @param kset A vector of integers specifying the factor/loadings pairs to be
#'   plotted. If \code{kset = NULL}, then all will be plotted.
#'
#' @param pm_which Whether to plot loadings \eqn{L} or factors \eqn{F} in the
#'   plots of posterior means.
#'   This parameter is ignored when \code{include_pm = FALSE}.
#'
#' @param pm_subset A vector of row indices \eqn{i} or column indices
#'   \eqn{j} (depending on the argument to \code{pm_which})
#'   specifying which values \eqn{\ell_{i \cdot}} or \eqn{f_{j \cdot}} are
#'   to be shown in the plots of posterior means. If the dataset has row or
#'   column names, then names rather than indices may be specified. If
#'   \code{pm_subset = NULL}, then all values will be plotted.
#'   This parameter is ignored when \code{include_pm = FALSE}.
#'
#' @param pm_groups A vector specifying the group to which each row of the data
#'   \eqn{y_{i \cdot}} or column \eqn{y_{\cdot j}} belongs
#'   (groups may be numeric indices or strings). If \code{pm_groups = NULL},
#'   then a bar plot of the ungrouped data is produced (see \code{include_pm}
#'   above). Otherwise, a group must be provided for each plotted row \eqn{i} or
#'   column \eqn{j}, so that
#'   the length of \code{pm_groups} is exactly equal to the number of rows or
#'   columns in the full dataset or, if \code{pm_subset} is specified, in the
#'   subsetted dataset. When \code{pm_groups} is not \code{NULL}, a set of
#'   overlapping histograms is produced for each factor/loadings pair, with
#'   one histogram per group (again see \code{include_pm}).
#'   This parameter is ignored when \code{include_pm = FALSE}.
#'
#' @param pm_colors A vector specifying a color for each bar (if
#'   \code{pm_groups = NULL}) or histogram (if \code{pm_groups} is not
#'   \code{NULL}). Passed directly to parameter \code{values} in \strong{ggplot2}
#'   function \code{\link[ggplot2]{scale_color_manual}}.
#'   This parameter is ignored when \code{include_pm = FALSE}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @return If arguments \code{include_scree} and \code{include_pm} specify that
#'   only one figure be produced, then \code{plot.flash()} returns a
#'   \code{ggplot2} object. If both figures are to be produced, then
#'   \code{plot.flash()} prints both plots but does not return a value.
#'
#' @method plot flash
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate left_join
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_col
#' @importFrom ggplot2 geom_histogram geom_vline after_stat
#' @importFrom ggplot2 scale_x_continuous scale_y_log10
#' @importFrom ggplot2 scale_fill_manual scale_color_manual
#' @importFrom ggplot2 scale_fill_brewer scale_color_brewer
#' @importFrom ggplot2 labs guides facet_wrap
#' @importFrom ggplot2 theme theme_minimal element_blank
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom tidyr pivot_longer
#' @importFrom stats density
#' @importFrom grDevices devAskNewPage
#'
#' @export
#'
plot.flash <- function(x,
                       include_scree = TRUE,
                       include_pm = TRUE,
                       order_by_pve = TRUE,
                       kset = NULL,
                       pm_which = c("factors", "loadings"),
                       pm_subset = NULL,
                       pm_groups = NULL,
                       pm_colors = NULL,
                       ...) {
  pm_which <- match.arg(pm_which)

  if (x$n_factors == 0) {
    stop("Flash object has no factors, so there is nothing to plot.")
  }
  if (!include_scree && !include_pm) {
    stop("Include either the scree plot (via argument include_scree) or the plot ",
         " of posterior means (via include_pm).")
  }

  # Bind variables to get rid of annoying R CMD check note:
  pve <- k.order <- Name <- k <- grp <- NULL

  # Scree plot:
  all.plots <- list()
  if (is.null(kset)) {
    kset <- 1:x$n_factors
  } else {
    must.be.valid.kset(get.fit(x), kset)
  }

  pve.df <- data.frame(
    k = kset,
    pve = x$pve[kset],
    k.order = kset
  )
  if (order_by_pve) {
    pve.df <- pve.df %>%
      mutate(k.order = rank(-pve))
  }

  if (include_scree) {
    # Include ~4 x-axis ticks with the distance a multiple of 5:
    K <- length(kset)
    if (K < 5) {
      tick_dist <- 1
    } else if (K < 10) {
      tick_dist <- 2
    } else {
      tick_dist <- max(1, floor(K / 20)) * 5
    }
    p1 <- ggplot(pve.df) +
      geom_line(aes(x = k.order, y = pve), color = "grey") +
      geom_point(aes(x = k.order, y = pve), color = "dodgerblue") +
      scale_x_continuous(breaks = seq(tick_dist, K, by = tick_dist)) +
      scale_y_log10() +
      labs(title = "Proportion of variance explained per factor") +
      labs(x = "k", y = "PVE") +
      theme_minimal()
    all.plots[["scree"]] <- p1
  }

  if (include_pm) {
    if (pm_which == "factors") {
      which.dim <- "column"
      val <- ldf(x, type = "i")$F
      colnames(val) <- 1:x$n_factors
    } else {
      which.dim <- "row"
      val <- ldf(x, type = "i")$L
      colnames(val) <- 1:x$n_factors
    }
    if (is.null(pm_subset)) {
      val <- val[, kset, drop = FALSE]
    } else {
      if (!all(pm_subset %in% 1:nrow(val)) &&
          !all(pm_subset %in% rownames(val))) {
        stop("Argument to pm_subset must be a vector of valid ", which.dim,
             " indices or valid ", which.dim, " names.")
      }
      val <- val[pm_subset, kset, drop = FALSE]
    }

    pm.df <- data.frame(val = val) %>%
      rownames_to_column(var = "Name")
    if (length(kset) > 1) {
      pm.df <- pm.df %>%
        pivot_longer(-Name, names_to = "k", values_to = "val") %>%
        mutate(k = as.numeric(str_remove(k, "val.")))
    } else {
      pm.df <- pm.df %>%
        mutate(k = kset)
    }
    pm.df <- pm.df %>%
      left_join(pve.df, by = "k")

    if (is.null(pm_groups)) {
      if (nrow(val) > 100) {
        warning("Consider setting argument pm_groups to produce a more",
                " readable plot.")
      }
      if (is.null(pm_colors)) {
        p2 <- ggplot(pm.df, aes(x = Name, y = val)) +
          geom_col(fill = "dodgerblue")
      } else {
        if (length(pm_colors) < nrow(val)) {
          stop("Argument to pm_colors must be a vector consisting of one ",
               "color for each ", which.dim, " in the data (or each ",
               "subsetted ", which.dim, ").")
        }
        p2 <- ggplot(pm.df, aes(x = Name, y = val, fill = Name)) +
          geom_col() +
          scale_fill_manual(values = pm_colors) +
          guides(fill = "none")
      }
      p2 <- p2 +
        facet_wrap(~k.order) +
        theme_minimal() +
        theme(axis.text = element_blank()) +
        labs(title = paste0("Posterior means for ", pm_which)) +
        labs(x = "", y = "")
    } else {
      if (length(pm_groups) != nrow(val)) {
        stop("Argument to pm_groups must be a vector with length equal to ",
             "the number of ", which.dim, "s in the data (or the number of ",
             "subsetted ", which.dim, "s.")
      }
      pm.df <- pm.df %>%
        mutate(grp = rep(factor(pm_groups), each = length(kset)))
      p2 <- ggplot(pm.df,
                   aes(x = val, y = after_stat(density), color = grp, fill = grp)) +
        geom_histogram(position = "identity", bins = 20, alpha = 0.5)
      if (is.null(pm_colors)) {
        if (length(unique(pm_groups)) > 9) {
          warning("Consider reducing the number of groups or defining a",
                  " custom color palette to produce a more readable plot.")
        } else {
          p2 <- p2 +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")
        }
      } else {
        if (length(pm_colors) < length(unique(pm_groups))) {
          stop("Argument to pm_colors must be a vector consisting of one ",
               "color for each unique group in argument pm_groups.")
        }
        p2 <- p2 +
          scale_color_manual(values = pm_colors) +
          scale_fill_manual(values = pm_colors)
      }
      p2 <- p2 +
        facet_wrap(~k.order, scales = "free_y") +
        geom_vline(xintercept = 0, color = "darkgrey") +
        theme_minimal() +
        theme(axis.text = element_blank()) +
        guides(color = "none") +
        labs(title = paste0("Posterior means for ", pm_which)) +
        labs(x = "", y = "", fill = "Group")
    }

    all.plots[["pm"]] <- p2
  }

  if (length(all.plots) == 1) {
    return(all.plots[[1]])
  } else {
    oask <- devAskNewPage(TRUE)
    print(all.plots[["scree"]])
    print(all.plots[["pm"]])
    devAskNewPage(oask)
    return(invisible())
  }
}

#' Fitted method for flash objects
#'
#' Given a \code{\link{flash}} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @return The matrix of "fitted values."
#'
#' @importFrom stats fitted
#' @method fitted flash
#'
#' @export
#'
fitted.flash <- function(object, ...) {
  return(fitted.flash_fit(get.fit(object)))
}

#' Fitted method for flash fit objects
#'
#' Given a \code{\link{flash_fit}} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash_fit}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @return The matrix of "fitted values."
#'
#' @importFrom stats fitted
#' @method fitted flash_fit
#'
#' @export
#'
fitted.flash_fit <- function(object, ...) {
  if (get.n.factors(object) == 0) {
    warning("Flash object does not have any factors.")
    return(matrix(data = 0, nrow = nrow(object$Y), ncol = ncol(object$Y)))
  }
  if (get.dim(object) > 2) {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }

  return(do.call(tcrossprod, get.EF(object)))
}

#' Residuals method for flash objects
#'
#' Given a \code{\link{flash}} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash
#'
#' @return The matrix of expected residuals.
#'
#' @importFrom stats residuals
#' @method residuals flash
#'
#' @export
#'
residuals.flash <- function(object, ...) {
  return(residuals.flash_fit(get.fit(object)))
}

#' Residuals method for flash fit objects
#'
#' Given a \code{\link{flash_fit}} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash_fit
#'
#' @return The matrix of expected residuals.
#'
#' @importFrom stats residuals
#' @method residuals flash_fit
#'
#' @export
#'
residuals.flash_fit <- function(object, ...) {
  if (uses.R(object)) {
    R <- get.R(object)
  } else {
    R <- get.Y(object) - lowrank.expand(get.EF(object))
  }

  if (any_missing(object)) {
    R[get.nonmissing(object) == 0] <- NA
  }

  return(R)
}

#' LDF method for flash and flash fit objects
#'
#' Given a \code{\link{flash}} or \code{\link{flash_fit}} object, returns the LDF
#'   decomposition \eqn{Y \approx LDF'}.
#'
#' When the prior families \eqn{G_\ell^{(k)}} and \eqn{G_f^{(k)}} are closed
#'   under scaling (as is typically the case), then the EBMF model (as
#'   described in the documention to function \code{\link{flash}}) is only
#'   identifiable up to scaling by a diagonal matrix \eqn{D}:
#'   \deqn{Y = LDF' + E.}
#'
#' Method \code{ldf} scales columns \eqn{\ell_k} and \eqn{f_k}
#'   so that, depending on the argument to parameter \code{type}, their
#'   1-norms, 2-norms, or infinity norms are equal to 1.
#'
#' @param object An object inheriting from class \code{flash} or
#'   \code{flash_fit}.
#'
#' @param type Takes identical arguments to function \code{\link[base]{norm}}. Use
#'   \code{"f"} or \code{"2"} for the 2-norm (Euclidean norm); \code{"o"} or
#'   \code{"1"} for the 1-norm (taxicab norm); and \code{"i"} or \code{"m"} for
#'   the infinity norm (maximum norm).
#'
#' @return A list with fields \code{L}, \code{D}, and \code{F}, each of which
#'   corresponds to one of the matrices in the decomposition \eqn{Y \approx LDF'}
#'   (with the columns of \eqn{L} and \eqn{F} scaled according to
#'   argument \code{type}). Note that \code{D} is returned as a vector rather
#'   than a matrix (the vector of diagonal entries in \eqn{D}). Thus, "fitted
#'   values" \eqn{LDF'} can be recovered as \code{L \%*\% diag(D) \%*\% t(F)}.
#'
#' @export
#'
ldf <- function(object, type) {
  UseMethod("ldf", object)
}

#' @describeIn ldf LDF decomposition for \code{\link{flash}} objects
#'
#' @method ldf flash
#'
#' @export
#'
ldf.flash <- function(object, type = "f") {
  return(ldf.flash_fit(get.fit(object), type = type))
}

#' @describeIn ldf LDF decomposition for \code{\link{flash_fit}} objects
#'
#' @method ldf flash_fit
#'
#' @export
#'
ldf.flash_fit <- function(object, type = "f") {
  type <- tolower(type)
  if (type == "2") {
    type <- "f"
  } else if (type == "1") {
    type <- "o"
  } else if (type == "m") {
    type <- "i"
  }

  if (get.n.factors(object) == 0) {
    stop("Flash fit does not have any factors.")
  }
  if (get.dim(object) > 2) {
    stop("S3 method \"ldf\" not available for tensors.")
  }

  ldf <- calc.normalized.loadings(object, type = type)

  ret <- list()
  ret$L <- ldf$normalized.loadings[[1]]
  ret$D <- ldf$scale.constants
  ret$F <- ldf$normalized.loadings[[2]]

  return(ret)
}

calc.normalized.loadings <- function(flash, for.pve = FALSE, type = "f") {
  ret <- list()

  if (for.pve) {
    loadings <- get.EF2(flash)
    norms <- lapply(loadings, colSums)
  } else {
    loadings <- get.EF(flash)
    if (type == "f") {
      norms <- lapply(loadings, function(x) {sqrt(colSums(x^2))})
    } else if (type == "o") {
      norms <- lapply(loadings, function(x) {colSums(abs(x))})
    } else if (type == "i") {
      norms <- lapply(loadings, function(x) {apply(abs(x), 2, max)})
    } else {
      stop("Norm type not recognized.")
    }
  }

  # Zero factors are "normalized" to zero.
  norms <- lapply(norms, function(x) {x[is.zero(flash)] <- Inf; x})
  L <- mapply(loadings, norms, FUN = function(X, y) {
    X / rep(y, each = nrow(X))
  }, SIMPLIFY = FALSE)
  # if (!for.pve) {
  #   L2 <- mapply(get.EF2(flash), norms, FUN = function(X, y) {
  #     X / rep(y^2, each = nrow(X))
  #   }, SIMPLIFY = FALSE)
  #   SD <- mapply(L2, L,
  #                FUN = function(EX2, EX) {sqrt(pmax(EX2 - EX^2, 0))},
  #                SIMPLIFY = FALSE)
  # }

  L <- propagate.names(L, flash)

  norms <- do.call(rbind, norms)
  ret$scale.constants <- apply(norms, 2, prod)
  ret$scale.constants[is.zero(flash)] <- 0
  ret$normalized.loadings <- L
  # if (!for.pve)
  #   ret$loading.SDs <- SD

  return(ret)
}

propagate.names <- function(mats, flash) {
  data.dimnames <- get.dimnames(flash)
  for (n in 1:get.dim(flash)) {
    if (!is.null(mats[[n]]) && !is.null(data.dimnames) && !is.null(data.dimnames[[n]]))
      rownames(mats[[n]]) <- data.dimnames[[n]]
  }
  return(mats)
}

#' @export
print.flash = function(x, ...) {
  if (x$n_factors == 0) {
    cat("Flash object with zero factors.\n")
  } else if (x$n_factors == 1) {
    cat("Flash object with one factor.\n")
    cat(sprintf("  Proportion of variance explained: %0.3f\n", x$pve))
  } else {
    cat(sprintf("Flash object with %d factors.\n", x$n_factors))
    if (all(x$pve < 0.001)) {
      cat("  All factors have PVE < 0.001.")
    } else {
      cat("  Proportion of variance explained")
      if (any(x$pve < 0.001)) {
        cat("*")
      }
      cat(":\n")
      for (k in 1:length(x$pve)) {
        if (x$pve[k] > 0.001) {
          cat(sprintf("    Factor %d: %0.3f\n", k, x$pve[k]))
        }
      }
      if (any(x$pve < 0.001)) {
        cat(paste("    *Factors with PVE < 0.001 are omitted from this",
                  "summary.\n"))
      }
    }
  }

  if (!is.na(x$elbo)) {
    cat(sprintf("  Variational lower bound: %0.3f\n", x$elbo))
  }
}
