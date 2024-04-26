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
                       order_by_pve = FALSE,
                       kset = NULL,
                       pm_which = c("factors", "loadings"),
                       pm_subset = NULL,
                       pm_groups = NULL,
                       pm_colors = NULL,
                       plot_type = c("scree",
                                     "bar",
                                     "heatmap",
                                     "histogram",
                                     "structure",
                                     "volcano",
                                     "data.frame"),
                       ...) {
  # TODO: use lifecycle package?
  if (!missing(include_scree) || !missing(include_pm)) {
    if (!missing(plot_type)) {
      warning("Please note that parameters 'include_scree' and 'include_pm' ",
              "have been soft-deprecated and will be removed in a future version ",
              "of flashier. Since 'plot_type' has been specified, 'include_scree' ",
              "and 'include_pm' will be ignored.")
    } else {
      warning("Please note that parameters 'include_scree' and 'include_pm' ",
              "have been soft-deprecated and will be removed in a future version ",
              "of flashier. To change the type of plot produced, please specify ",
              "argument 'plot_type'.")
      if (include_pm & !include(include_scree)) {
        if (is.null(pm_groups)) {
          plot_type <- "bar"
        } else {
          plot_type <- "histogram"
        }
      }
    }
  }

  if (plot_type == "scree") {
    ret <- flash_plot_scree(
      x, order_by_pve, kset
    )
  } else if (plot_type == "bar") {
    ret <- flash_plot_bar(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors
    )
  } else if (plot_type == "heatmap") {
    ret <- flash_plot_heatmap(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "histogram") {
    ret <- flash_plot_histogram(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "structure") {
    ret <- flash_plot_structure(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "volcano") {
    ret <- flash_plot_volcano(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "data.frame") {
    ret <- flash_plot_dataframe(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors
    )
  }

  return(ret)
}

#' @importFrom ggplot2 ggplot aes geom_line geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_log10 labs
#' @importFrom dplyr group_by summarize
#' @importFrom cowplot theme_cowplot
#'
flash_plot_scree <- function(fl,
                             order_by_pve = FALSE,
                             kset = NULL) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = "factors",
                             pm_subset = NULL,
                             pm_groups = NULL,
                             pm_colors = NULL)
  df <- df |>
    group_by(k_order, pve) |>
    summarize(.groups = "drop")

  # Bind variables to get rid of annoying R CMD check note:
  pve <- k_order <- NULL

  p <- ggplot(df) +
    geom_line(aes(x = k_order, y = pve), color = "grey") +
    geom_point(aes(x = k_order, y = pve), color = "dodgerblue") +
    scale_y_log10() +
    labs(title = "Scree plot") +
    labs(x = "k", y = "PVE") +
    theme_cowplot(font_size = 10)

  # Include ~4 x-axis ticks with the distance a multiple of 5:
  K <- length(unique(df$k_order))
  if (K < 5) {
    tick_dist <- 1
  } else if (K < 10) {
    tick_dist <- 2
  } else {
    tick_dist <- max(1, floor(K / 20)) * 5
  }

  p <- p +
    scale_x_continuous(breaks = seq(tick_dist, K, by = tick_dist))

  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_col
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggplot2 facet_wrap element_blank
#' @importFrom cowplot theme_cowplot
#'
flash_plot_bar <- function(fl,
                           order_by_pve = FALSE,
                           kset = NULL,
                           pm_which = c("factors", "loadings"),
                           pm_subset = NULL,
                           pm_groups = NULL,
                           pm_colors = NULL) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = pm_colors)

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
    theme_cowplot(font_size = 10) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "")
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_histogram after_stat
#' @importFrom ggplot2 scale_fill_brewer scale_color_brewer
#' @importFrom ggplot2 scale_fill_identity scale_color_identity
#' @importFrom ggplot2 facet_wrap geom_vline
#' @importFrom ggplot2 theme element_blank
#' @importFrom ggplot2 guides labs
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#'
flash_plot_histogram <- function(fl,
                                 order_by_pve = FALSE,
                                 kset = NULL,
                                 pm_which = c("factors", "loadings"),
                                 pm_subset = NULL,
                                 pm_groups = NULL,
                                 pm_colors = NULL,
                                 bins = 20,
                                 alpha = 0.5) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = if (is.null(pm_groups)) NULL else pm_colors)

  # Bind variables to get rid of annoying R CMD check note:
  val <- group <- color <- k_order <- NULL

  if (is.null(df$group)) {
    p <- ggplot(df, aes(x = val, y = after_stat(density))) +
      geom_histogram(position = "identity", bins = bins, fill = "dodgerblue")
  } else {
    p <- ggplot(df, aes(x = val, y = after_stat(density),
                        color = color, fill = color)) +
      geom_histogram(position = "identity", bins = bins, alpha = alpha) +
      scale_color_identity() +
      scale_fill_identity()
  }

  p <- p +
    facet_wrap(~k_order, scales = "free_y") +
    geom_vline(xintercept = 0, color = "darkgrey") +
    theme_cowplot(font_size = 10) +
    theme(axis.text = element_blank()) +
    guides(color = "none") +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "", fill = "Group")
  return(p)
}

#' @importFrom ggplot2 labs
#' @importFrom fastTopics structure_plot
#'
flash_plot_structure <- function(fl,
                                 order_by_pve = FALSE,
                                 kset = NULL,
                                 pm_which = c("factors", "loadings"),
                                 pm_subset = NULL,
                                 pm_groups = NULL,
                                 pm_colors = NULL,
                                 gap = 1,
                                 ...) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = NULL)

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
                        gap = gap,
                        ...)
  } else {
    p <- structure_plot(Lmat,
                        topics = rev(unique(df$k_order)),
                        grouping = group,
                        colors = pm_colors,
                        gap = gap,
                        ...)
  }

  if (any(df$val < 0) && any(df$val > 0)) {
    warning("Structure plots were designed to visualize sets of nonnegative ",
            "memberships or loadings. Structure plots that include negative ",
            "loadings are often difficult to interpret, so a heatmap should ",
            "typically be preferred when visualizing a combination of negative ",
            "and positive loadings.")
  }

  p <- p +
    labs(y = "loading", color = "factor", fill = "factor")
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_y_continuous labs
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#'
flash_plot_heatmap <- function(fl,
                               order_by_pve = FALSE,
                               kset = NULL,
                               pm_which = c("factors", "loadings"),
                               pm_subset = NULL,
                               pm_groups = NULL,
                               pm_colors = c("darkblue", "white", "darkred"),
                               gap = 1,
                               ...) {
  # Use flash_plot_structure to get embedding:
  struct_p <- flash_plot_structure(fl = fl,
                                   order_by_pve = order_by_pve,
                                   kset = kset,
                                   pm_which = pm_which,
                                   pm_subset = pm_subset,
                                   pm_groups = pm_groups,
                                   pm_colors = NULL,
                                   gap = gap,
                                   ...)
  struct_df <- struct_p$data

  # Retrieve group information:
  struct_ticks <- struct_p$plot_env$ticks

  # Topics get reversed by plot_structure; re-reverse them:
  struct_df$topic <- factor(struct_df$topic, level = rev(levels(struct_df$topic)))

  p <- ggplot(struct_df, aes(x = topic, y = sample, fill = prop)) +
    geom_tile(width = 0.8) +
    scale_fill_gradient2(low = pm_colors[3], mid = pm_colors[2], high = pm_colors[1]) +
    scale_y_continuous(breaks = struct_ticks, labels = names(struct_ticks)) +
    labs(x = "factor", y = "", fill = "loading") +
    theme_cowplot(font_size = 10)
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom ggplot2 scale_color_identity
#' @importFrom ggplot2 facet_wrap labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by mutate summarize arrange left_join
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#'
flash_plot_volcano <- function(fl,
                               order_by_pve = FALSE,
                               kset = NULL,
                               pm_which = c("factors", "loadings"),
                               pm_subset = NULL,
                               pm_groups = NULL,
                               pm_colors = NULL,
                               n_labels = 0,
                               label_size = 3,
                               ...) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = pm_colors)

  fl <- flash_fit(fl)
  Y <- get.Y(fl)
  if (match.arg(pm_which) == "loadings") {
    n <- 1
    other_n <- 2
    which_dim <- "row"
  } else {
    n <- 2
    other_n <- 1
    which_dim <- "column"
  }
  dimsums <- nmode.prod.r1(Y, r1.ones(fl), n)
  dimmeans <- dimsums / get.data.dims(Y)[other_n]
  mean_df <- data.frame(
    name = get.data.dimnames(Y)[[n]],
    mean_exp = dimmeans
  )
  df <- df |>
    left_join(mean_df, by = "name")

  if (is.null(df$color)) {
    df$color = "black"
  }

  if (n_labels > 0) {
    df <- df |>
      group_by(k_order) |>
      mutate(label = ifelse(rank(-abs(val)) > n_labels, "", name),
             color = ifelse(label != "", "dodgerblue", "gray80"))
  }

  p <- ggplot(df, aes(x = val, y = mean_exp, color = color)) +
    geom_point() +
    facet_wrap(~k_order) +
    labs(x = "loading", y = "size factor") +
    theme_cowplot(font_size = 10)

  if (is.null(pm_groups)) {
    p <- p +
      scale_color_identity()
  } else {
    color_df <- df |>
      group_by(group, color) |>
      summarize(.groups = "drop") |>
      arrange(group)
    p <- p +
      scale_color_identity(guide = "legend",
                           name = "",
                           labels = color_df$group,
                           breaks = color_df$color)
  }

  if (n_labels > 0) {
    p <- p +
      geom_text_repel(aes(label = label),
                      size = label_size,
                      ...)
  }

  return(p)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom Polychrome kelly.colors
#'
flash_plot_dataframe <- function(fl,
                                 order_by_pve = FALSE,
                                 kset = NULL,
                                 pm_which = c("factors", "loadings"),
                                 pm_subset = NULL,
                                 pm_groups = NULL,
                                 pm_colors = NULL) {
  if (fl$n_factors == 0) {
    stop("Flash object has no factors, so there is nothing to plot.")
  }

  pm_which <- match.arg(pm_which)

  if (is.null(kset)) {
    kset <- 1:fl$n_factors
  } else {
    must.be.valid.kset(get.fit(fl), kset)
  }

  pve <- fl$pve[kset]
  if (order_by_pve) {
    k_order <- rank(-pve)
  } else {
    k_order <- 1:length(kset)
  }

  if (pm_which == "factors") {
    which.dim <- "column"
    val <- ldf(fl, type = "i")$F
    colnames(val) <- 1:fl$n_factors
  } else {
    which.dim <- "row"
    val <- ldf(fl, type = "i")$L
    colnames(val) <- 1:fl$n_factors
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
    if (!is.null(pm_colors)) {
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
    n_colors <- length(unique(pm_groups))
    if (!is.null(pm_colors)) {
      if (length(pm_colors) != n_colors || !is.character(pm_colors)) {
        stop("When 'pm_groups' is set, 'pm_colors' must be a character vector ",
             "consisting of one color for each unique value in 'pm_groups'.")
      }
    } else {
      if (n_colors < 10) {
        pm_colors <- brewer.pal(n_colors, "Set1")
      } else {
        pm_colors <- rep(kelly.colors(22)[-1], length.out = n_colors)
      }
    }
    df$color <- rep(pm_colors[factor(pm_groups)], ncol(val))
  }

  return(df)
}
