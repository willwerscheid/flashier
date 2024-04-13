#' Plot method for flash objects
#'
#' Given a \code{\link{flash}} object, produces up to two figures: one showing
#'   the proportion of variance explained per factor/loadings pair, and one that
#'   plots posterior means for either factors or loadings (depending on the
#'   argument to parameter \code{pm_which}).
#'
#' @param x An object inheriting from class \code{flash}.
#'
#' @param include_scree This parameter has been soft-deprecated; please use
#'   \code{plot_type} instead.
#'
#' @param include_pm This parameter has been soft-deprecated; please use
#'   \code{plot_type} instead.
#'
#'   Whether to include a figure showing the posterior means for
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
#' @param pm_which Whether to plot loadings \eqn{L} or factors \eqn{F}. This
#'   parameter is ignored when \code{plot_type = "scree"}.
#'
#' @param pm_subset A vector of row indices \eqn{i} or column indices
#'   \eqn{j} (depending on the argument to \code{pm_which})
#'   specifying which values \eqn{\ell_{i \cdot}} or \eqn{f_{j \cdot}} are
#'   to be shown. If the dataset has row or column names, then names rather
#'   than indices may be specified. If \code{pm_subset = NULL}, then all values
#'   will be plotted. This parameter is ignored when \code{plot_type = "scree"}.
#'
#' @param pm_groups A vector specifying the group to which each row of the data
#'   \eqn{y_{i \cdot}} or column \eqn{y_{\cdot j}} belongs (groups may be
#'   numeric indices or strings). A group must be provided for each plotted row
#'   \eqn{i} or column \eqn{j}, so that the length of \code{pm_groups} is
#'   exactly equal to the number of rows or columns in the full dataset or, if
#'   \code{pm_subset} is specified, in the subsetted dataset. For effects, see
#'   parameter \code{plot_type}.
#'
#' @param pm_colors A character vector specifying a color for each unique group
#'   specified by \code{pm_groups}, or, if \code{pm_groups = NULL}, a vector
#'   specifying a color for each plotted row \eqn{i} or column \eqn{j}. For
#'   effects, see parameter \code{plot_type}.
#'
#' @param plot_type The type of plot to return. Options include:
#'   \describe{
#'     \item{\code{"scree"}}{A scree plot showing the proportion of variance
#'       explained per factor/loadings pair.}
#'     \item{\code{"bar"}}{A simple bar plot of posterior means for loadings or
#'       factors (depending on argument \code{pm_which}), with one bar per
#'       row or column. Colors of bars are specified by argument
#'       \code{pm_colors}. This type of plot is most useful when rows or columns
#'       are small in number or ordered in a logical fashion (e.g., spatially).}
#'     \item{\code{"heatmap"}}{A heatmap showing posterior means for loadings or
#'       factors, with rows or columns grouped according to \code{pm_groups} and
#'       arranged within groups using the embedding provided by \code{fastTopics}
#'       function \code{\link[fastTopics]{structure_plot}}. Note that
#'       heatmaps ignore argument \code{pm_colors}.}
#'     \item{\code{"histogram"}}{Overlapping semi-transparent histograms of
#'       posterior means for loadings or factors, with one histogram per group
#'       specified by \code{pm_groups} (or a single histogram if \code{pm_groups}
#'       is \code{NULL}). Colors of histograms are specified by \code{pm_colors}.}
#'     \item{\code{"structure"}}{A "structure plot" (stacked bar plot) produced
#'       using function \code{\link[fastTopics]{structure_plot}} in package
#'       \code{fastTopics}. Here \code{pm_colors} specifies the colors of
#'       different factor/loadings pairs (as specified by \code{kset}) rather
#'       than different groups (as specified by \code{pm_groups}). Note that
#'       factors/loadings must be nonnegative for structure plots to make
#'       sense.}
#'   }
#'
#' @param ... Additional parameters to be passed to function
#'   \code{\link[fastTopics]{structure_plot}} when \code{plot_type} is
#'   \code{"heatmap"} or \code{"structure"}. In particular, the \code{gap}
#'   parameter is useful for establishing visual separation among groups.
#'
#' @return If arguments \code{include_scree} and \code{include_pm} specify that
#'   only one figure be produced, then \code{plot.flash()} returns a
#'   \code{ggplot2} object. If both figures are to be produced, then
#'   \code{plot.flash()} prints both plots but does not return a value.
#'
#' @method plot flash
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
                       plot_type = c("scree",
                                     "bar",
                                     "heatmap",
                                     "histogram",
                                     "structure"),
                       ...) {
  pm_which <- match.arg(pm_which)
  plot_type <- match.arg(plot_type)

  if (x$n_factors == 0) {
    stop("Flash object has no factors, so there is nothing to plot.")
  }

  if (is.null(kset)) {
    kset <- 1:x$n_factors
  } else {
    must.be.valid.kset(get.fit(x), kset)
  }

  pve <- x$pve[kset]
  if (order_by_pve) {
    k_order <- rank(-pve)
  } else {
    k_order <- kset
  }

  if (plot_type == "scree") {
    df <- data.frame(
      k = kset,
      pve = pve,
      k_order = k_order
    )
  } else {
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
      pm_subset <- 1:nrow(val)
    } else if (!all(pm_subset %in% 1:nrow(val)) &&
               !all(pm_subset %in% rownames(val))) {
        stop("Argument to pm_subset must be a vector of valid ", which.dim,
             " indices or valid ", which.dim, " names.")
    }
    val <- val[pm_subset, kset, drop = FALSE]

    df <- data.frame(
      name = rep(rownames(val), ncol(val)),
      val = as.vector(val),
      k = rep(kset, each = nrow(val)),
      pve = rep(pve, each = nrow(val)),
      k_order = rep(k_order, each = nrow(val))
    )
    if (is.null(pm_groups)) {
      if (!is.null(pm_colors) && plot_type != "structure") {
        if (length(pm_colors) != nrow(val) || !is.character(pm_colors)) {
          stop("When 'pm_groups' is NULL, 'pm_colors' must a character vector ",
               "consisting of one color for each ", which.dim, " in the data ",
               "(or each subsetted ", which.dim, ").")
        }
        df$color <- rep(pm_colors, ncol(val))
      }
    } else {
      if (length(pm_groups) != nrow(val)) {
        stop("'pm_groups' must be NULL or a vector consisting of one value ",
             "(numeric or character) for each ", which.dim, " in the data ",
             "(or each subsetted ", which.dim, ").")
      }
      df$group <- rep(pm_groups, ncol(val))
      if (!is.null(pm_colors) && plot_type != "structure") {
        if (length(pm_colors) != length(unique(pm_groups))
            || !is.character(pm_colors)) {
          stop("When 'pm_groups' is set, 'pm_colors' must be a character vector ",
               "consisting of one color for each unique value in 'pm_groups'.")
        }
      }
      df$color <- rep(pm_colors[factor(pm_groups)], ncol(val))
    }
  }

  # TODO: use lifecycle package?
  if (!missing(include_scree) || !missing(include_pm)) {
    warning("Please note that parameters 'include_scree' and 'include_pm' ",
            "have been soft-deprecated and will be removed in a future version ",
            "of flashier. To change the type of plot produced, please specify ",
            "argument 'plot_type'.")
  }

  if (plot_type == "scree") {
    p <- plot_scree(df)
  } else if (plot_type == "bar") {
    p <- plot_bar(df, pm_which)
  } else if (plot_type == "heatmap") {
    p <- plot_heatmap(df, pm_which, ...)
  } else if (plot_type == "histogram") {
    p <- plot_histogram(df, pm_which)
  } else if (plot_type == "structure") {
    p <- plot_structure(df, pm_which, pm_colors, ...)
  }

  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_line geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_log10
#' @importFrom ggplot2 labs theme_minimal
#'
plot_scree <- function(df) {
  # Bind variables to get rid of annoying R CMD check note:
  pve <- k_order <- NULL

  # Include ~4 x-axis ticks with the distance a multiple of 5:
  K <- nrow(df)
  if (K < 5) {
    tick_dist <- 1
  } else if (K < 10) {
    tick_dist <- 2
  } else {
    tick_dist <- max(1, floor(K / 20)) * 5
  }
  p <- ggplot(df) +
    geom_line(aes(x = k_order, y = pve), color = "grey") +
    geom_point(aes(x = k_order, y = pve), color = "dodgerblue") +
    scale_x_continuous(breaks = seq(tick_dist, K, by = tick_dist)) +
    scale_y_log10() +
    labs(title = "Proportion of variance explained per factor") +
    labs(x = "k", y = "PVE") +
    theme_minimal()
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_col
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggplot2 facet_wrap theme_minimal element_blank
#'
plot_bar <- function(df, pm_which) {
  # Bind variables to get rid of annoying R CMD check note:
  name <- val <- color <- k_order <- NULL

  if (is.null(df$color)) {
    p <- ggplot(df, aes(x = name, y = val)) +
      geom_col(fill = "dodgerblue")
  } else {
    p <- ggplot(df, aes(x = name, y = val, fill = color)) +
      geom_col() +
      scale_fill_identity()
  }
  p <- p +
    facet_wrap(~k_order) +
    theme_minimal() +
    theme(axis.text = element_blank()) +
    labs(title = paste0("Posterior means for ", pm_which)) +
    labs(x = "", y = "")
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_y_continuous labs
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#'
plot_heatmap <- function(df, pm_which, ...) {
  # Use plot_structure to get embedding:
  struct_p <- plot_structure(df, pm_which, rep("black", nrow(df)), ...)
  struct_df <- struct_p$data

  # Retrieve group information:
  struct_ticks <- struct_p$plot_env$ticks

  # Topics get reversed by plot_structure; re-reverse them:
  struct_df$topic <- factor(struct_df$topic, level = rev(levels(struct_df$topic)))

  p <- ggplot(struct_df, aes(x = topic, y = sample, fill = prop)) +
    geom_tile(width = 0.8) +
    scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue") +
    scale_y_continuous(breaks = struct_ticks, labels = names(struct_ticks)) +
    labs(x = "factor", y = "", fill = "membership") +
    theme_cowplot(font_size = 10)
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_histogram after_stat
#' @importFrom ggplot2 scale_fill_brewer scale_color_brewer
#' @importFrom ggplot2 scale_fill_identity scale_color_identity
#' @importFrom ggplot2 facet_wrap geom_vline
#' @importFrom ggplot2 theme theme_minimal element_blank
#' @importFrom ggplot2 guides labs
#' @importFrom stats density
#'
plot_histogram <- function(df, pm_which) {
  # Bind variables to get rid of annoying R CMD check note:
  val <- group <- color <- k_order <- NULL

  if (is.null(df$group)) {
    p <- ggplot(df, aes(x = val, y = after_stat(density))) +
      geom_histogram(position = "identity", bins = 20, fill = "dodgerblue")
  } else {
    p <- ggplot(df, aes(x = val, y = after_stat(density),
                        color = color, fill = color)) +
      geom_histogram(position = "identity", bins = 20, alpha = 0.5)

    if (is.null(df$color)) {
      if (length(unique(df$group)) > 9) {
        warning("Consider reducing the number of groups or defining a ",
                " custom color palette via 'pm_colors' to produce a more readable ",
                "plot.")
      } else {
        p <- p +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1")
      }
    } else {
      p <- p +
        scale_color_identity() +
        scale_fill_identity()
    }
  }

  p <- p +
    facet_wrap(~k_order, scales = "free_y") +
    geom_vline(xintercept = 0, color = "darkgrey") +
    theme_minimal() +
    theme(axis.text = element_blank()) +
    guides(color = "none") +
    labs(title = paste0("Posterior means for ", pm_which)) +
    labs(x = "", y = "", fill = "Group")
  return(p)
}

#' @importFrom ggplot2 labs
#' @importFrom fastTopics structure_plot
#'
plot_structure <- function(df, pm_which, pm_colors, ...) {
  Lmat <- matrix(df$val, ncol = length(unique(df$k)))
  colnames(Lmat) <- paste0("k", unique(df$k))
  if (!is.null(df$group)) {
    group <- df$group[1:nrow(Lmat)]
  } else {
    group <- rep("", nrow(Lmat))
  }
  if (is.null(pm_colors)) {
    p <- structure_plot(Lmat,
                        topics = rev(unique(df$k_order)),
                        grouping = group,
                        ...)
  } else {
    p <- structure_plot(Lmat,
                        topics = rev(unique(df$k_order)),
                        grouping = group,
                        colors = pm_colors,
                        ...)
  }
  p <- p +
    labs(y = "membership", color = "factor", fill = "factor")
  return(p)
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
