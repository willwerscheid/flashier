#' Plot method for flash objects
#'
#' Plots a \code{\link{flash}} object. Several types of plot are possible:
#'   see parameter \code{plot_type} below as well as functions
#'   \code{\link{flash_plot_scree}}, \code{\link{flash_plot_bar}},
#'   \code{\link{flash_plot_heatmap}}, \code{\link{flash_plot_histogram}},
#'   \code{\link{flash_plot_scatter}}, and \code{\link{flash_plot_structure}}.
#'
#' @param x An object inheriting from class \code{flash}.
#'
#' @param include_scree `r lifecycle::badge("deprecated")` This parameter has been deprecated; please use
#'   \code{plot_type} instead.
#'
#' @param include_pm `r lifecycle::badge("deprecated")` This parameter has been deprecated; please use
#'   \code{plot_type} instead.
#'
#' @param order_by_pve If \code{order_by_pve = TRUE}, then factor/loadings pairs
#'   will be ordered according to proportion of variance explained, from
#'   highest to lowest. (By default, they are plotted in the same order as
#'   \code{kset}; or, if \code{kset} is \code{NULL}, then they are plotted in
#'   the same order as they are found in \code{fl}.)
#'
#' @param kset A vector of integers specifying the factor/loadings pairs to be
#'   plotted. If \code{order_by_pve = FALSE}, then \code{kset} also specifies the
#'   \emph{order} in which they are to be plotted.
#'
#' @param pm_which Whether to plot loadings \eqn{L} or factors \eqn{F}.
#'
#' @param pm_subset A vector of row indices \eqn{i} or column indices
#'   \eqn{j} (depending on the argument to \code{pm_which})
#'   specifying which values \eqn{\ell_{i \cdot}}{l_{i.}} or \eqn{f_{j \cdot}}{f_{j.}} are
#'   to be shown. If the dataset has row or column names, then names rather
#'   than indices may be specified. If \code{pm_subset = NULL}, then all values
#'   will be plotted.
#'
#' @param pm_groups A vector specifying the group to which each row of the data
#'   \eqn{y_{i \cdot}}{y_{i.}} or column \eqn{y_{\cdot j}}{y_{.j}} belongs (groups may be
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
#'       explained per factor/loadings pair. See
#'       \code{\link{flash_plot_scree}}.}
#'     \item{\code{"bar"}}{A bar plot of posterior means for loadings or
#'       factors (depending on argument \code{pm_which}), with one bar per
#'       row or column. Colors of bars are specified by argument
#'       \code{pm_colors}. This type of plot is most useful when rows or columns
#'       are small in number or ordered in a logical fashion (e.g., spatially).
#'       See \code{\link{flash_plot_bar}}.}
#'     \item{\code{"heatmap"}}{A heatmap showing posterior means for loadings or
#'       factors, with rows or columns grouped using a 1-d embedding. Here
#'       \code{pm_color} specifies the diverging color gradient (low-mid-high).
#'       See \code{\link{flash_plot_heatmap}}.}
#'     \item{\code{"histogram"}}{Overlapping semi-transparent histograms of
#'       posterior means for loadings or factors, with one histogram per group
#'       specified by \code{pm_groups} (or a single histogram if \code{pm_groups}
#'       is \code{NULL}). Colors of histograms are specified by \code{pm_colors}.
#'       See \code{\link{flash_plot_histogram}}.}
#'     \item{\code{"scatter"}}{A scatter plot showing the relationship between
#'       posterior means for loadings or factors and a user-supplied covariate.
#'       If a covariate is not supplied, then data column or row means will be
#'       used. Colors of points are specified by \code{pm_colors}. See
#'       \code{\link{flash_plot_scatter}}.}
#'     \item{\code{"structure"}}{A "structure plot" (stacked bar plot) produced
#'       using function \code{\link[fastTopics]{structure_plot}} in package
#'       \code{fastTopics}. Here \code{pm_colors} specifies the colors of
#'       different factor/loadings pairs (as specified by \code{kset}) rather
#'       than different groups (as specified by \code{pm_groups}). Note that
#'       factors/loadings must be nonnegative for structure plots to make
#'       sense. See \code{\link{flash_plot_structure}}.}
#'   }
#'
#' @param ... Additional parameters to be passed to respective
#'   \code{flash_plot_xxx} functions. See
#'   \code{\link{flash_plot_scree}}, \code{\link{flash_plot_bar}},
#'   \code{\link{flash_plot_heatmap}}, \code{\link{flash_plot_histogram}},
#'   \code{\link{flash_plot_scatter}}, and \code{\link{flash_plot_structure}}
#'   for details.
#'
#' @return A \code{ggplot} object.
#'
#' @method plot flash
#'
#' @importFrom lifecycle deprecate_warn
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
                                     "scatter",
                                     "structure"),
                       ...) {
  deprecation_details <- paste(
    "Parameters 'include_scree' and 'include_pm' will be removed in a future",
    "release. To change the type of plot produced, please specify 'plot_type'",
    "instead."
  )
  if (!missing(include_scree)) {
    deprecate_warn(
      when = "1.0.28",
      what = "plot.flash(include_scree)",
      details = deprecation_details
    )
  }
  if (!missing(include_pm)) {
    deprecate_warn(
      when = "1.0.28",
      what = "plot.flash(include_scree)",
      details = deprecation_details
    )
  }

  if (!missing(include_scree) || !missing(include_pm)) {
    if (!missing(plot_type)) {
      plot_type <- match.arg(plot_type)
    } else {
      if (include_pm && (missing(include_scree) || !include_scree)) {
        if (is.null(pm_groups)) {
          plot_type <- "bar"
        } else {
          plot_type <- "histogram"
        }
      } else {
        plot_type <- "scree"
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
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "heatmap") {
    ret <- flash_plot_heatmap(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "histogram") {
    ret <- flash_plot_histogram(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "scatter") {
    ret <- flash_plot_scatter(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
    )
  } else if (plot_type == "structure") {
    ret <- flash_plot_structure(
      x, order_by_pve, kset, pm_which, pm_subset, pm_groups, pm_colors, ...
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
#' @param fl An object inheriting from class \code{flash}.
#'
#' @param labels Whether to label the points in the scree plot with the
#'   indices of the factor/loading pairs they correspond to. Labels appear
#'   as "k1", "k2", "k3", etc.
#'
#' @param label_size The size of the label text (in millimeters).
#'
#' @param max_overlaps A (nonnegative) integer. For each text label, the number
#'   of overlaps with other text labels or other data points are counted, and
#'   the text label is excluded if it has too many overlaps.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(gtex)
#' fl <- flash(gtex, greedy_Kmax = 4L, backfit = FALSE)
#' flash_plot_scree(fl)
#'
#' # For the full range of labelling options provided by the ggrepel package, set
#' #   labels = FALSE (the default setting) and add geom_text_repel() manually:
#' library(ggrepel)
#' flash_plot_scree(fl) + geom_text_repel(min.segment.length = 0)
#'
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_line geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_log10 labs
#' @importFrom cowplot theme_cowplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#'
#' @export
#'
flash_plot_scree <- function(fl,
                             order_by_pve = FALSE,
                             kset = NULL,
                             labels = FALSE,
                             label_size = 3,
                             max_overlaps = Inf) {
  df <- flash_plot_dataframe(fl = fl,
                             order_by_pve = order_by_pve,
                             kset = kset,
                             pm_which = "factors",
                             pm_subset = NULL,
                             pm_groups = NULL,
                             pm_colors = NULL)
  df <- df |>
    group_by(.data$pve, .data$k) |>
    summarize(.groups = "drop") |>
    mutate(k_numeric = as.numeric(.data$k))

  p <- ggplot(df, aes(x = .data$k_numeric, y = .data$pve, label = .data$k)) +
    geom_line(color = "grey") +
    geom_point(color = "dodgerblue") +
    scale_y_log10() +
    labs(x = "k", y = "PVE") +
    theme_cowplot(font_size = 10)

  # Include ~4 x-axis ticks with the distance a multiple of 5:
  K <- length(levels(df$k))
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
      geom_text_repel(size = label_size,
                      max.overlaps = max_overlaps) +
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
#' @inheritParams flash_plot_scree
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
#' @return A \code{ggplot} object.
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
#' @importFrom ggplot2 ggplot aes geom_col
#' @importFrom ggplot2 scale_fill_identity labs facet_wrap
#' @importFrom ggplot2 theme element_text element_blank
#' @importFrom cowplot theme_cowplot
#' @importFrom rlang .data
#'
#' @export
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

  if (is.null(df$color)) {
    p <- ggplot(df, aes(x = .data$name, y = .data$val)) +
      geom_col(fill = "dodgerblue")
  } else {
    p <- ggplot(df, aes(x = .data$name, y = .data$val, fill = .data$color)) +
      geom_col() +
      scale_fill_identity()
  }
  p <- p +
    theme_cowplot(font_size = 10) +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "")

  if (length(levels(df$k)) > 1) {
    p <- p +
      facet_wrap(~.data$k, ...)
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
#' @param ... Additional arguments to be passed to
#'   \code{\link[ggplot2]{facet_wrap}} (e.g., \code{nrow} or \code{ncol}).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom stats density
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat
#' @importFrom ggplot2 scale_fill_identity scale_color_identity
#' @importFrom ggplot2 facet_wrap geom_vline
#' @importFrom ggplot2 theme element_blank
#' @importFrom ggplot2 guides labs
#' @importFrom cowplot theme_cowplot
#' @importFrom rlang .data
#'
#' @export
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

  if (is.null(pm_groups)) {
    p <- ggplot(df, aes(x = .data$val, y = after_stat(density))) +
      geom_histogram(position = "identity", bins = bins, fill = "dodgerblue")
  } else {
    color_df <- df |>
      group_by(.data$group, .data$color) |>
      summarize(.groups = "drop") |>
      arrange(.data$group)

    p <- ggplot(df, aes(x = .data$val, y = after_stat(density),
                        fill = .data$color)) +
      geom_histogram(position = "identity", bins = bins, color = "white") +
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
    # geom_vline(xintercept = 0, color = "darkgrey") +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_blank()) +
    labs(title = paste0("Posterior means (", pm_which, ")")) +
    labs(x = "", y = "", fill = "Group")

  if (length(levels(df$k)) > 1) {
    p <- p +
      facet_wrap(~.data$k, scales = "free_y", ...)
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
#' @inheritParams flash_plot_scree
#'
#' @inheritParams flash_plot_bar
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
#' @param labels Whether to label the points with the largest (absolute)
#'   posterior means. If \code{labels = TRUE}, then \code{n_labels} points will
#'   be labelled using \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param n_labels A (nonnegative) integer. The number of points to label. If
#'   \code{n_labels} is set to a positive integer but \code{labels = FALSE},
#'   then the \code{n_labels} points with the largest (absolute) posterior
#'   means will be highlighted in blue but not labelled. This can be useful for
#'   tweaking labels using the full range of options provided by
#'   \code{\link[ggrepel]{geom_text_repel}}. For an example, see below.
#'
#' @param ... Additional arguments to be passed to
#'   \code{\link[ggplot2]{facet_wrap}} (e.g., \code{nrow} or \code{ncol}).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(gtex)
#' fl <- flash(gtex, greedy_Kmax = 4L, backfit = FALSE)
#' flash_plot_scatter(fl)
#'
#' # Label axes and points:
#' library(ggplot2)
#' flash_plot_scatter(fl, labels = TRUE, n_labels = 3) +
#'   labs(y = "mean z-score across all SNPs")
#'
#' # For the full range of labelling options provided by the ggrepel package, set
#' #   labels = FALSE (the default setting) and add geom_text_repel() manually:
#' library(ggrepel)
#' flash_plot_scatter(fl, labels = FALSE, n_labels = 3) +
#'   geom_text_repel(size = 2.5, min.segment.length = 0)
#'
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom ggplot2 scale_color_identity
#' @importFrom ggplot2 facet_wrap labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by mutate summarize arrange left_join
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @export
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
                               labels = FALSE,
                               n_labels = 0,
                               label_size = 3,
                               max_overlaps = Inf,
                               ...) {
  if (labels && missing(n_labels)) {
    n_labels <- 10
    message("The number of labels has been set to 10. To change this setting ",
            "and suppress this message, please set argument 'n_labels'.")
  }

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

  if (is.null(df$color)) {
    df$color = "black"
  }

  if (n_labels > 0) {
    df <- df |>
      group_by(.data$k) |>
      mutate(label = ifelse(rank(-abs(.data$val)) > n_labels, "", .data$name)) |>
      mutate(color = ifelse(.data$label != "", "dodgerblue", "gray80"))
  } else {
    df$label = ""
  }

  p <- ggplot(df, aes(x = .data$val, y = .data$covariate,
                      color = .data$color, label = .data$label)) +
    geom_point(shape = shape) +
    labs(x = paste(sub("s", "", pm_which), "posterior mean"),
         y = ylab) +
    theme_cowplot(font_size = 10)

  if (length(levels(df$k)) > 1) {
    p <- p +
      facet_wrap(~k, ...)
  }

  if (is.null(pm_groups)) {
    p <- p +
      scale_color_identity()
  } else {
    color_df <- df |>
      group_by(.data$group, .data$color) |>
      summarize(.groups = "drop") |>
      arrange(.data$group)
    p <- p +
      scale_color_identity(guide = "legend",
                           name = "",
                           labels = color_df$group,
                           breaks = color_df$color)
  }

  if (labels) {
    p <- p +
      geom_text_repel(size = label_size,
                      max.overlaps = max_overlaps)
  }

  return(p)
}

#' Create structure plot of factors or loadings for a flash fit
#'
#' Creates a \dQuote{structure plot} (stacked bar plot) of posterior means for factors
#'   \eqn{f_{jk}} or loadings \eqn{\ell_{ik}}{l_{ik}}. Different \dQuote{topics} or components
#'   (that is, the different factor/loadings pairs, as specified by \code{kset})
#'   are represented by different colors. Values are normalized so that the
#'   maximum absolute value for each factor \eqn{f_{\cdot k}}{f_{.k}} or set of
#'   loadings \eqn{\ell_{\cdot k}}{l_{.k}} is equal to 1 and then stacked (see
#'   \code{\link{ldf.flash}}). Note that structure plots were designed for
#'   nonnegative loadings or \dQuote{memberships}; if posterior means are not
#'   nonnegative then a different type of plot should be used (e.g.,
#'   \code{\link{flash_plot_heatmap}}). By default, a 1-d embedding is used to
#'   arrange the rows \eqn{i} or columns \eqn{j}. This step is usually essential
#'   to creating a readable structure plot; for details, see
#'   \code{\link[fastTopics]{structure_plot}}.
#'
#' @inheritParams plot.flash
#'
#' @inheritParams flash_plot_scree
#'
#' @param pm_colors The colors of the \dQuote{topics} or components (factor/loadings
#'   pairs).
#'
#' @param gap The horizontal spacing between groups. Ignored if \code{pm_groups}
#'   is not provided.
#'
#' @param ... Additional parameters to be passed to
#'   \code{\link[fastTopics]{structure_plot}}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom fastTopics structure_plot
#' @importFrom ggplot2 labs
#'
#' @export
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

  Lmat <- matrix(df$val, ncol = length(levels(df$k)))
  colnames(Lmat) <- levels(df$k)
  if (!is.null(df$group)) {
    group <- df$group[1:nrow(Lmat)]
  } else {
    group <- rep("", nrow(Lmat))
  }

  if (is.null(pm_colors)) {
    p <- structure_plot(Lmat,
                        topics = rev(levels(df$k)),
                        grouping = group,
                        gap = gap,
                        ...)
  } else {
    p <- structure_plot(Lmat,
                        topics = rev(levels(df$k)),
                        grouping = group,
                        colors = pm_colors,
                        gap = gap,
                        ...)
  }

  if (any(df$val < 0) && any(df$val > 0)) {
    warning("Structure plots were designed to visualize sets of nonnegative ",
            "memberships or loadings. Structure plots that include negative ",
            "values are often difficult to interpret, so a heatmap should ",
            "typically be preferred when visualizing a combination of negative ",
            "and positive values.")
  }

  p <- p +
    labs(y = "posterior mean", color = "component", fill = "component")
  return(p)
}

#' Create heatmap of factors or loadings for a flash fit
#'
#' Creates a heatmap of posterior means for factors \eqn{f_{jk}} or loadings
#'   \eqn{\ell_{ik}}{l_{ik}}. Values are normalized so that the maximum absolute value
#'   for each factor \eqn{f_{\cdot k}}{f_{.k}} or set of loadings \eqn{\ell_{\cdot k}}{l_{.k}}
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
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_y_continuous labs
#' @importFrom cowplot theme_cowplot
#' @importFrom stats density
#' @importFrom rlang .data
#'
#' @export
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
  struct_ticks <- struct_p$plot_env$ticks + 0.5

  # Topics get reversed by plot_structure; re-reverse them:
  struct_df$topic <- factor(struct_df$topic,
                            levels = rev(levels(struct_df$topic)))

  p <- ggplot(struct_df, aes(x = .data$topic, y = .data$sample, fill = .data$prop)) +
    geom_tile(width = 0.8, height = 1) +
    scale_fill_gradient2(low = pm_colors[3], mid = pm_colors[2], high = pm_colors[1]) +
    scale_y_continuous(breaks = struct_ticks, labels = names(struct_ticks)) +
    labs(x = "component", y = "", fill = "posterior\n mean") +
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
    k_order <- order(-pve)
  } else {
    k_order <- 1:length(kset)
  }
  k <- factor(paste0("k", kset),
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
    k = rep(k, each = nrow(val)),
    val = as.vector(val),
    pve = rep(pve, each = nrow(val))
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
