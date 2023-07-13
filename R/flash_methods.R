#' Plot method for flash objects
#'
#' Given a \code{\link{flash}} object, produces up to two figures: one showing
#'   the proportion of variance explained per factor, and one that plots
#'   posterior means.
#'
#' @param x An object inheriting from class \code{flash}.
#'
#' @param incl.scree Whether to include a figure showing the proportion of
#'   variance explained by each factor/loading pair ("scree plot").
#'
#' @param incl.pm Whether to include a figure showing the posterior means for
#'   either loadings \eqn{L} or factors \eqn{F} (depending on the argument to
#'   \code{pm.which}). If argument \code{pm.groups}
#'   is left unspecified, then bar plots will be produced, with each bar
#'   corresponding to a single value \eqn{\ell_{ik}} or \eqn{f_{jk}}.
#'   Otherwise, overlapping histograms will be
#'   produced, with each histogram corresponding to one of the groups
#'   specified by \code{pm.groups}. One plot is produced for each
#'   factor/loading pair \eqn{k}.
#'
#' @param order.by.pve If \code{TRUE}, then the factor/loading pairs will be
#'   re-ordered according to proportion of variance explained (from
#'   highest to lowest).
#'
#' @param kset A vector of integers specifying the factor/loading pairs to be
#'   plotted. If \code{kset = NULL}, then all will be plotted.
#'
#' @param pm.which Whether to plot loadings \eqn{L} or factors \eqn{F} in the
#'   plots of posterior means.
#'
#' @param pm.subset A vector of loading indices \eqn{i} or factor indices
#'   \eqn{j} (depending on the argument to \code{pm.which})
#'   specifying which values \eqn{\ell_{i \cdot}} or \eqn{f_{j \cdot}} are
#'   to be shown in the plots of posterior means. If the dataset has row or
#'   column names, then names rather than indices may be specified. If
#'   \code{pm.subset = NULL}, then all values will be plotted.
#'
#' @param pm.groups A vector specifying the group to which each row
#'   \eqn{y_{i \cdot}} or column \eqn{y_{\cdot j}} of the data belongs
#'   (groups may be numeric indices or strings). If \code{pm.groups = NULL},
#'   then a bar plot of the ungrouped data will be produced (see \code{incl.pm}
#'   above). Otherwise, a group must be provided for each observation, so that
#'   the length of \code{pm.groups} is exactly equal to the number of rows or
#'   columns in the full dataset or, if \code{pm.subset} is specified, in the
#'   subsetted dataset. If \code{pm.groups} is not \code{NULL}, then a set of
#'   overlapping histograms will be produced for each factor/loading pair, with
#'   one histogram per group (again see \code{incl.pm}).
#'
#' @param pm.colors A vector specifying a color for each bar (if
#'   \code{pm.groups = NULL}) or histogram (if \code{pm.groups} is not
#'   \code{NULL}). Passed directly to argument \code{values} in \strong{ggplot2}
#'   function \code{\link[ggplot2]{scale_color_manual}}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @seealso \code{\link{flash}}
#'
#' @method plot flash
#'
#' @importFrom magrittr `%>%`
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
#' @importFrom grDevices devAskNewPage
#'
#' @export
#'
plot.flash <- function(x,
                       incl.scree = TRUE,
                       incl.pm = TRUE,
                       order.by.pve = TRUE,
                       kset = NULL,
                       pm.which = c("factors", "loadings"),
                       pm.subset = NULL,
                       pm.groups = NULL,
                       pm.colors = NULL,
                       ...) {
  pm.which <- match.arg(pm.which)

  if (x$n.factors == 0) {
    stop("Flash object has no factors, so there is nothing to plot.")
  }
  if (!incl.scree && !incl.pm) {
    stop("Include either the scree plot (via argument incl.scree) or the plot ",
         " of posterior means (via incl.pm).")
  }

  # Scree plot:
  all.plots <- list()
  if (is.null(kset)) {
    kset <- 1:x$n.factors
  } else {
    must.be.valid.kset(get.fit(x), kset)
  }

  pve.df <- data.frame(
    k = kset,
    pve = x$pve[kset],
    k.order = kset
  )
  if (order.by.pve) {
    pve.df <- pve.df %>%
      mutate(k.order = rank(-pve))
  }

  if (incl.scree) {
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

  if (incl.pm) {
    if (pm.which == "factors") {
      which.dim <- "column"
      val <- ldf(x, type = "i")$F
      colnames(val) <- 1:x$n.factors
    } else {
      which.dim <- "row"
      val <- ldf(x, type = "i")$L
      colnames(val) <- 1:x$n.factors
    }
    if (is.null(pm.subset)) {
      val <- val[, kset, drop = FALSE]
    } else {
      if (!all(pm.subset %in% 1:nrow(val)) &&
          !all(pm.subset %in% rownames(val))) {
        stop("Argument to pm.subset must be a vector of valid ", which.dim,
             " indices or valid ", which.dim, " names.")
      }
      val <- val[pm.subset, kset, drop = FALSE]
    }

    pm.df <- data.frame(val = val) %>%
      rownames_to_column(var = "Name") %>%
      pivot_longer(-Name, names_to = "k", values_to = "val") %>%
      mutate(k = as.numeric(str_remove(k, "val."))) %>%
      left_join(pve.df, by = "k")

    if (is.null(pm.groups)) {
      if (nrow(val) > 100) {
        warning("Consider setting argument pm.groups to produce a more",
                " readable plot.")
      }
      if (is.null(pm.colors)) {
        p2 <- ggplot(pm.df) +
          geom_col(aes(x = Name, y = val), fill = "dodgerblue")
      } else {
        if (length(pm.colors) < nrow(val)) {
          stop("Argument to pm.colors must be a vector consisting of one ",
               "color for each ", which.dim, " in the data (or each ",
               "subsetted ", which.dim, ").")
        }
        p2 <- ggplot(pm.df) +
          geom_col(aes(x = Name, y = val, fill = Name)) +
          scale_fill_manual(values = pm.colors) +
          guides(fill = "none")
      }
      p2 <- p2 +
        facet_wrap(~k.order) +
        theme_minimal() +
        theme(axis.text = element_blank()) +
        labs(title = paste0("Posterior means for ", pm.which)) +
        labs(x = "", y = "")
    } else {
      if (length(pm.groups) != nrow(val)) {
        stop("Argument to pm.groups must be a vector with length equal to ",
             "the number of ", which.dim, "s in the data (or the number of ",
             "subsetted ", which.dim, "s.")
      }
      pm.df <- pm.df %>%
        mutate(grp = rep(factor(pm.groups), each = length(kset)))
      p2 <- ggplot(pm.df) +
        geom_histogram(
          aes(x = val, y = after_stat(density), color = grp, fill = grp),
          position = "identity", bins = 20, alpha = 0.5
        )
      if (is.null(pm.colors)) {
        if (length(unique(pm.groups)) > 9) {
          warning("Consider reducing the number of groups or defining a",
                  " custom color palette to produce a more readable plot.")
        } else {
          p2 <- p2 +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")
        }
      } else {
        if (length(pm.colors) < length(unique(pm.groups))) {
          stop("Argument to pm.colors must be a vector consisting of one ",
               "color for each unique group in argument pm.groups.")
        }
        p2 <- p2 +
          scale_color_manual(values = pm.colors) +
          scale_fill_manual(values = pm.colors)
      }
      p2 <- p2 +
        facet_wrap(~k.order, scales = "free_y") +
        geom_vline(xintercept = 0, color = "darkgrey") +
        theme_minimal() +
        theme(axis.text = element_blank()) +
        guides(color = "none") +
        labs(title = paste0("Posterior means for ", pm.which)) +
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
#' Given a \code{flash} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @seealso \code{\link{flash}}
#'
#' @method fitted flash
#'
#' @export
#'
fitted.flash <- function(object, ...) {
  return(fitted.flash.fit(get.fit(object)))
}

#' Fitted method for flash fit objects
#'
#' Given a \code{flash.fit} object, returns the "fitted values"
#'   \eqn{E(LF') = E(L) E(F)'}.
#'
#' @param object An object inheriting from class \code{flash.fit}.
#'
#' @param ... Additional parameters are ignored.
#'
#' @seealso \code{\link{flash}}
#'
#' @method fitted flash.fit
#'
#' @export
#'
fitted.flash.fit <- function(object, ...) {
  if (get.n.factors(object) == 0) {
    stop("Flash object does not have any factors.")
  }
  if (get.dim(object) > 2) {
    stop("S3 method \"fitted\" not yet implemented for tensors.")
  }

  return(do.call(tcrossprod, get.EF(object)))
}

#' Residuals method for flash objects
#'
#' Given a \code{flash} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash
#'
#' @seealso \code{\link{flash}}
#'
#' @method residuals flash
#'
#' @export
#'
residuals.flash <- function(object, ...) {
  return(residuals.flash.fit(get.fit(object)))
}

#' Residuals method for flash fit objects
#'
#' Given a \code{flash.fit} object, returns the expected residuals
#'   \eqn{Y - E(LF') = Y - E(L) E(F)'}.
#'
#' @inheritParams fitted.flash.fit
#'
#' @seealso \code{\link{flash}}
#'
#' @method residuals flash.fit
#'
#' @export
#'
residuals.flash.fit <- function(object, ...) {
  if (uses.R(object)) {
    R <- get.R(object)
  } else {
    R <- get.Y(object) - lowrank.expand(get.EF(object))
  }

  if (any.missing(object)) {
    R[get.nonmissing(object) == 0] <- NA
  }

  return(R)
}

#' LDF method for flash and flash fit objects
#'
#' Given a \code{flash} or \code{flash.fit} object, returns the LDF
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
#'   \code{flash.fit}.
#'
#' @param type Takes identical arguments to function \code{\link[base]{norm}}. Use
#'   \code{"f"} or \code{"2"} for the 2-norm (Euclidean norm); \code{"o"} or
#'   \code{"1"} for the 1-norm (taxicab norm); and \code{"i"} or \code{"m"} for
#'   the infinity norm (maximum norm).
#'
#' @seealso \code{\link{flash}}
#'
#' @export
#'
ldf <- function(object, type) {
  UseMethod("ldf", object)
}

#' @describeIn ldf LDF decomposition for flash objects
#'
#' @method ldf flash
#'
#' @export
#'
ldf.flash <- function(object, type = "f") {
  return(ldf.flash.fit(get.fit(object), type = type))
}

#' @describeIn ldf LDF decomposition for flash fit objects
#'
#' @method ldf flash.fit
#'
#' @export
#'
ldf.flash.fit <- function(object, type = "f") {
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
  if (x$n.factors == 0) {
    cat("Flash object with zero factors.\n")
  } else if (x$n.factors == 1) {
    cat("Flash object with one factor.\n")
    cat(sprintf("  Proportion of variance explained: %0.3f\n", x$pve))
  } else {
    cat(sprintf("Flash object with %d factors.\n", x$n.factors))
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
