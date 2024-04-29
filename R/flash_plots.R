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
#' @param order_by_pve If \code{TRUE}, then the factor/loadings pairs will be
#'   re-ordered according to proportion of variance explained (from
#'   highest to lowest).
#'
#' @param kset A vector of integers specifying the factor/loadings pairs to be
#'   plotted. If \code{kset = NULL}, then all will be plotted.
#'
#' @param pm_which Whether to plot loadings \eqn{L} or factors \eqn{F}.
#'
#' @param pm_subset A vector of row indices \eqn{i} or column indices
#'   \eqn{j} (depending on the argument to \code{pm_which})
#'   specifying which values \eqn{\ell_{i \cdot}} or \eqn{f_{j \cdot}} are
#'   to be shown. If the dataset has row or column names, then names rather
#'   than indices may be specified. If \code{pm_subset = NULL}, then all values
#'   will be plotted.
#'
#' @param pm_groups A vector specifying the group to which each row of the data
#'   \eqn{y_{i \cdot}} or column \eqn{y_{\cdot j}} belongs (groups may be
#'   numeric indices or strings). A group must be provided for each plotted row
#'   \eqn{i} or column \eqn{j}, so that the length of \code{pm_groups} is
#'   exactly equal to the number of rows or columns in the full dataset or, if
#'   \code{pm_subset} is specified, in the subsetted dataset.
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
#' @param ... Additional parameters to be passed to particular
#'   \code{flash_plot_xxx} functions. See individual functions for details.
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
                                     "scatter",
                                     "data.frame"),
                       ...) {
  # TODO: use lifecycle package?
  if (!missing(include_scree) || !missing(include_pm)) {
    if (!missing(plot_type)) {
      warning("Please note that parameters 'include_scree' and 'include_pm' ",
              "have been soft-deprecated and will be removed in a future version ",
              "of flashier. Since 'plot_type' has been specified, 'include_scree' ",
              "and 'include_pm' will be ignored.")
      plot_type <- match.arg(plot_type)
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
  } else {
    plot_type <- match.arg(plot_type)
  }

  if (plot_type == "scree") {
    ret <- flash_plot_scree(
      x, order_by_pve, kset, ...
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
  } else if (plot_type == "scatter") {
    ret <- flash_plot_scatter(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "data.frame") {
    ret <- flash_plot_dataframe(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors
    )
  }

  return(ret)
}

#' Create a scree plot for a flash fit
#'
#' A scree plot is a line plot showing the proportion of variance explained by
#'   each factor/loadings pair in a flash fit. Note that since EBMF does not
#'   require factors and loadings to be orthogonal, "PVE" should be
#'   interpreted loosely: for example, the total proportion of variance
#'   explained could be larger than 1.
#'
#' Unlike scree plots for PCA, a scree plot for a flash fit is not in general
#'   monotonically decreasing. To ensure a monotonically decreasing scree
#'   plot, set \code{order_by_pve = TRUE}. Note, however, that if this is done
#'   then the numbers on the \eqn{x}-axis will no longer match the indices of
#'   the components in the flash fit. This can also be true if argument
#'   \code{kset} has been specified. Thus one should consider setting
#'   \code{labels = TRUE} when \code{order_by_pve = TRUE} or when \code{kset}
#'   has been specified.
#'
#' @inheritParams plot.flash
#'
#' @param labels Whether to label the points in the scree plot with the
#'   indices of the factor/loading pairs they correspond to. Labels appear
#'   as "k1", "k2", "k3", etc.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_line geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_log10 labs
#' @importFrom cowplot theme_cowplot
#' @importFrom ggrepel geom_text_repel
#'
flash_plot_scree <- function(fl,
                             order_by_pve = FALSE,
                             kset = NULL,
                             labels = FALSE) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = "factors",
                             pm_subset = NULL,
                             pm_groups = NULL,
                             pm_colors = NULL)
  # Bind variables to get rid of annoying R CMD check note:
  k_order <- pve <- k_factor <- NULL

  df <- df |>
    group_by(k_order, pve, k_factor) |>
    summarize(.groups = "drop")

  p <- ggplot(df, aes(x = k_order, y = pve)) +
    geom_line(color = "grey") +
    geom_point(color = "dodgerblue") +
    scale_y_log10() +
    labs(title = "Scree plot", x = "k", y = "PVE") +
    theme_cowplot(font_size = 10)

  # Include ~4 x-axis ticks with the distance a multiple of 5:
  K <- length(levels(df$k_factor))
  if (K < 5) {
    tick_dist <- 1
  } else if (K < 10) {
    tick_dist <- 2
  } else {
    tick_dist <- max(1, floor(K / 20)) * 5
  }
  p <- p +
    scale_x_continuous(breaks = seq(tick_dist, K, by = tick_dist))

  if (labels) {
    p <- p +
      geom_text_repel(aes(label = k_factor)) +
      theme(axis.text.x = element_blank())
  }

  return(p)
}

#' Create bar plots of factors or loadings for a flash fit
#'
#' Creates a bar plot or sequence of bar plots, one for each value of \eqn{k} in
#'   \code{kset}, with bars corresponding to individual posterior means for factors
#'   \eqn{f_{jk}} or loadings \eqn{\ell_{ik}}. Values are normalized so that the
#'   maximum absolute value for each factor \eqn{f_{\cdot k}} or set of
#'   loadings \eqn{\ell_{\cdot k}} is equal to 1 (see \code{\link{ldf.flash}}).
#'   This type of plot is most useful when rows \eqn{i = 1, \ldots, n} or columns
#'   \eqn{j = 1, \ldots, p} are small in number or ordered in a logical fashion
#'   (e.g., spatially).
#'
#' When there is more than one value of \eqn{k} in \code{kset}, a sequence of
#'   panels is created using the \code{\link[ggplot2]{facet_wrap}} function from
#'   the \code{ggplot2} package. In each panel, the order of bars is determined
#'   by the order of the corresponding rows or columns in the data matrix;
#'   they can be re-arranged using the \code{pm_subset} argument.
#'
#' @inheritParams plot.flash
#'
#' @param pm_colors A character vector specifying a color for each unique group
#'   specified by \code{pm_groups}, or, if \code{pm_groups = NULL}, a vector
#'   specifying a color for each plotted row \eqn{i} or column \eqn{j}. Defines
#'   the color (fill) of the bars.
#'
#' @param labels Whether to label the bars along the \eqn{x}-axis. The
#'   appearance of the labels (size, angle, etc.) can be adjusted using
#'   \code{ggplot2}'s theme system; see below for an example.
#'
#' @param ... Additional arguments to be passed to
#'   \code{\link[ggplot2]{facet_wrap}} (e.g., \code{nrow} or \code{ncol}).
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_col
#' @importFrom ggplot2 scale_fill_identity labs facet_wrap
#' @importFrom ggplot2 theme element_text element_blank
#' @importFrom cowplot theme_cowplot
#'
#' @examples
#' data(gtex)
#' fl <- flash(gtex, greedy_Kmax = 4L, backfit = FALSE)
#' flash_plot_bar(fl, pm_colors = gtex_colors)
#'
#' # Tweaks are often required to get x-axis labels to look good:
#' library(ggplot2)
#' flash_plot_bar(fl, pm_colors = gtex_colors, labels = TRUE, ncol = 1) +
#'   theme(axis.text.x = element_text(size = 8, angle = 60))
#'
flash_plot_bar <- function(fl,
                           order_by_pve = FALSE,
                           kset = NULL,
                           pm_which = c("factors", "loadings"),
                           pm_subset = NULL,
                           pm_groups = NULL,
                           pm_colors = NULL,
                           labels = FALSE,
                           ...) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = pm_colors)

  # Bind variables to get rid of annoying R CMD check note:
  name <- val <- color <- k_factor <- NULL

  if (is.null(df$color)) {
    p <- ggplot(df, aes(x = name, y = val)) +
      geom_col(fill = "dodgerblue")
  } else {
    p <- ggplot(df, aes(x = name, y = val, fill = color)) +
      geom_col() +
      scale_fill_identity()
  }
  p <- p +
    theme_cowplot(font_size = 10) +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "")

  if (length(levels(df$k_factor)) > 1) {
    p <- p +
      facet_wrap(~k_factor, ...)
  }

  if (labels) {
    p <- p +
      theme(axis.text.x = element_text(
        angle = 80, hjust = 1, size = 8
      ))
  } else {
    p <- p +
      theme(axis.text.x = element_blank()) +
      theme(axis.ticks.x = element_blank())
  }

  return(p)
}

#' Create histograms of factors or loadings for a flash fit
#'
#' Creates a histogram or sequence of histograms of posterior means for factors
#'   \eqn{f_{jk}} or loadings \eqn{\ell_{ik}}. One plot is created for each
#'   value of \eqn{k} in \code{kset}. Values are normalized so that the
#'   maximum absolute value for each factor \eqn{f_{\cdot k}} or set of
#'   loadings \eqn{\ell_{\cdot k}} is equal to 1 (see \code{\link{ldf.flash}}).
#'   If \code{pm_groups} is specified, then overlapping semi-transparent
#'   histograms are created, with one histogram per group specified by
#'   \code{pm_groups}. This option works best when the number of groups is small
#'   or when groups are well separated across components.
#'
#' @inheritParams flash_plot_bar
#'
#' @param pm_colors A character vector specifying a color for each unique group
#'   specified by \code{pm_groups}. Defines the color and fill of the histograms.
#'
#' @param binwidth The width of the bins (a numeric value). The default is to
#'   use the number of bins in \code{bins}, covering the range of the data.
#'
#' @param bins Number of bins. Overriden by \code{binwidth}. Defaults to 30.
#'
#' @param alpha A transparency value between 0 (transparent) and 1 (opaque).
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom stats density
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat
#' @importFrom ggplot2 scale_fill_identity scale_color_identity
#' @importFrom ggplot2 facet_wrap geom_vline
#' @importFrom ggplot2 theme element_blank
#' @importFrom ggplot2 guides labs
#' @importFrom cowplot theme_cowplot
#'
flash_plot_histogram <- function(fl,
                                 order_by_pve = FALSE,
                                 kset = NULL,
                                 pm_which = c("factors", "loadings"),
                                 pm_subset = NULL,
                                 pm_groups = NULL,
                                 pm_colors = NULL,
                                 binwidth = NULL,
                                 bins = NULL,
                                 alpha = 0.5,
                                 ...) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = pm_which,
                             pm_subset = pm_subset,
                             pm_groups = pm_groups,
                             pm_colors = if (is.null(pm_groups)) NULL else pm_colors)

  # Bind variables to get rid of annoying R CMD check note:
  val <- group <- color <- k_order <- NULL

  if (is.null(pm_groups)) {
    p <- ggplot(df, aes(x = val, y = after_stat(density))) +
      geom_histogram(position = "identity", bins = bins, fill = "dodgerblue")
  } else {
    color_df <- df |>
      group_by(group, color) |>
      summarize(.groups = "drop") |>
      arrange(group)

    p <- ggplot(df, aes(x = val, y = after_stat(density),
                        color = color, fill = color)) +
      geom_histogram(position = "identity", bins = bins, alpha = alpha) +
      scale_color_identity(guide = "legend",
                           name = "",
                           labels = color_df$group,
                           breaks = color_df$color) +
      scale_fill_identity(guide = "legend",
                          name = "",
                          labels = color_df$group,
                          breaks = color_df$color)
  }

  p <- p +
    geom_vline(xintercept = 0, color = "darkgrey") +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_blank()) +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "", fill = "Group")

  if (length(levels(df$k_factor)) > 1) {
    p <- p +
      facet_wrap(~k_factor, scales = "free_y", ...)
  }
  return(p)
}

#' Create scatter plots of factors or loadings for a flash fit
#'
#' Creates a scatter plot or sequence of scatter plots, with position along the
#'   \eqn{x}-axis defined by posterior means for factors \eqn{f_{jk}} or loadings
#'   \eqn{\ell_{ik}} and position along the \eqn{y}-axis defined by a
#'   user-supplied covariate. If a covariate is not supplied, then plots will
#'   use data column or row means, \eqn{\frac{1}{n} \sum_{i = 1}^n y_{ij}} or
#'   \eqn{\frac{1}{p} \sum_{j = 1}^p y_{ij}}. One plot is created for
#'   each value of \eqn{k} in \code{kset}. Values are normalized so that the
#'   maximum absolute value for each factor \eqn{f_{\cdot k}} or set of
#'   loadings \eqn{\ell_{\cdot k}} is equal to 1 (see \code{\link{ldf.flash}}).
#'
#' @inheritParams plot.flash
#'
#' @param pm_colors A character vector specifying a color for each unique group
#'   specified by \code{pm_groups}, or, if \code{pm_groups = NULL}, a vector
#'   specifying a color for each plotted row \eqn{i} or column \eqn{j}. Defines
#'   the colors of the points.
#'
#' @param covariate A numeric vector with one value for each plotted row \eqn{i}
#'   or column \eqn{j}. These values are mapped onto the plots' \eqn{y}-axis.
#'
#' @param shape The symbol used for the plots' points. See
#'   \code{\link[ggplot2]{aes_linetype_size_shape}}.
#'
#' @param n_labels A (nonnegative) integer. If \code{n_labels > 0}, then the
#'   points with the \code{n_labels} largest (absolute) posterior means will be
#'   labelled using \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param label_size The size of the label text (in millimeters).
#'
#' @param ... Additional arguments to be passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom ggplot2 scale_color_identity
#' @importFrom ggplot2 facet_wrap labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by mutate summarize arrange left_join
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#'
flash_plot_scatter <- function(fl,
                               order_by_pve = FALSE,
                               kset = NULL,
                               pm_which = c("factors", "loadings"),
                               pm_subset = NULL,
                               pm_groups = NULL,
                               pm_colors = NULL,
                               covariate = NULL,
                               shape = 1,
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

  if (!is.null(covariate)) {
    df$covariate <- rep(covariate, length.out = nrow(df))
    ylab <- "covariate"
  } else {
    # Use row/column data means as the default covariate.
    if (match.arg(pm_which) == "loadings") {
      n <- 1
      other_n <- 2
      which_dim <- "row"
    } else {
      n <- 2
      other_n <- 1
      which_dim <- "column"
    }

    fl <- flash_fit(fl)
    Y <- get.Y(fl)
    dimsums <- nmode.prod.r1(Y, r1.ones(fl), n)
    dimmeans <- dimsums / get.data.dims(Y)[other_n]
    covariate_df <- data.frame(
      name = get.data.dimnames(Y)[[n]],
      covariate = dimmeans
    )
    df <- df |>
      left_join(covariate_df, by = "name")

    ylab <- paste("data", which_dim, "mean")
  }

  # Bind variables to get rid of annoying R CMD check note:
  val <- group <- color <- k_factor <- NULL

  if (is.null(df$color)) {
    df$color = "black"
  }

  if (n_labels > 0) {
    df <- df |>
      group_by(k_factor) |>
      mutate(label = ifelse(rank(-abs(val)) > n_labels, "", name),
             color = ifelse(label != "", "dodgerblue", "gray80"))
  }

  p <- ggplot(df, aes(x = val, y = covariate, color = color)) +
    geom_point(shape = shape) +
    labs(x = "loading", y = ylab) +
    theme_cowplot(font_size = 10)

  if (length(levels(df$k_factor)) > 1) {
    p <- p +
      facet_wrap(~k_factor)
  }

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

#' Create structure plot of factors or loadings for a flash fit
#'
#' Creates a "structure plot" (stacked bar plot) of posterior means for factors
#'   \eqn{f_{jk}} or loadings \eqn{\ell_{ik}}. Different "topics" or components
#'   (that is, the different factor/loadings pairs, as specified by \code{kset})
#'   are represented by different colors. Values are normalized so that the
#'   maximum absolute value for each factor \eqn{f_{\cdot k}} or set of
#'   loadings \eqn{\ell_{\cdot k}} is equal to 1 and then stacked (see
#'   \code{\link{ldf.flash}}). Note that structure plots were designed for
#'   nonnegative loadings or "memberships"; if posterior means are not
#'   nonnegative then a different type of plot should be used (e.g.,
#'   \code{\link{flash_plot_heatmap}}). By default, a 1-d embedding is used to
#'   arrange the rows \eqn{i} or columns \eqn{j}. This step is usually essential
#'   to creating a readable structure plot; for details, see
#'   \code{\link[fastTopic]{structure_plot}}.
#'
#' @inheritParams plot.flash
#'
#' @param pm_colors The colors of the "topics" or components (factor/loadings
#'   pairs).
#'
#' @param gap The horizontal spacing between groups. Ignored if \code{pm_groups}
#'   is not provided.
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link[fastTopics]{structure_plot}}.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom fastTopics structure_plot
#' @importFrom ggplot2 labs
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

#' Create heatmap of factors or loadings for a flash fit
#'
#' Creates a heatmap of posterior means for factors \eqn{f_{jk}} or loadings
#'   \eqn{\ell_{ik}}. Values are normalized so that the maximum absolute value
#'   for each factor \eqn{f_{\cdot k}} or set of loadings \eqn{\ell_{\cdot k}}
#'   is equal to 1 (see \code{\link{ldf.flash}}).
#'
#' By default, a 1-d embedding is used to arrange the rows \eqn{i} or columns
#'   \eqn{j} in a "smart" manner. This behavior can be overridden via argument
#'   \code{loadings_order}, which is passed to function
#'   \code{\link[fastTopics]{structure_plot}}.
#'
#' @inheritParams flash_plot_structure
#'
#' @param pm_colors A character vector of length 1, 2, or 3 defining the
#'   diverging color gradient (low-mid-high) to be used by the heatmap. The
#'   midpoint is set at zero. If one or two colors are supplied, then the
#'   "mid" color will be set to white. If one color is supplied, then the "low"
#'   and "high" colors (used for, respectively, negative and positive posterior
#'   means) will be the same. If two are supplied, then the "low" color should
#'   be provided first, followed by the "high" color. If all three are supplied,
#'   then the "low" color should be provided first, followed by the "mid" color,
#'   followed by the "high" color provided. The default color gradient is
#'   \code{darkblue} for "low" (negative posterior means), white for "mid"
#'   (zero), and \code{darkred} for "high" (positive posterior means).
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link[fastTopics]{structure_plot}} (which is primarily used to
#'   arrange the rows \eqn{i} or columns \eqn{j}).
#'
#' @return A \code{ggplot2} object.
#'
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
                               pm_colors = NULL,
                               gap = 1,
                               ...) {

  if (is.null(pm_colors)) {
    pm_colors <- c("darkred", "white", "darkblue")
  } else if (length(pm_colors) == 1) {
    pm_colors <- c(pm_colors, "white", pm_colors)
  } else if (length(pm_colors == 2)) {
    pm_colors <- c(pm_colors[1], "white", pm_colors[2])
  } else if (length(pm_colors) != 3) {
    stop("When creating a heatmap, pm_colors must be NULL or a vector of",
         "length 1, 2, or 3.")
  }

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
    geom_tile(width = 0.8, height = 1) +
    scale_fill_gradient2(low = pm_colors[3], mid = pm_colors[2], high = pm_colors[1]) +
    scale_y_continuous(breaks = struct_ticks, labels = names(struct_ticks)) +
    labs(x = "factor", y = "", fill = "loading") +
    theme_cowplot(font_size = 10)
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
  k_factor <- factor(paste0("k", kset),
                     levels = paste0("k", kset[k_order]))

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
    k_order = rep(k_order, each = nrow(val)),
    k_factor = rep(k_factor, each = nrow(val))
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
