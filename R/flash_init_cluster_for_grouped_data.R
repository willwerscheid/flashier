# TODO: handle fl in argument to init.
#' @importFrom tictoc tic toc
#' @importFrom snow makeCluster clusterExport clusterEvalQ sendCall
#' @importFrom snow checkForRemoteErrors recvResult
#' @export
flash_init_cluster_for_grouped_data <- function(dat,
                                                groups,
                                                nCores,
                                                K = 50L,
                                                ebnm_fn = ebnm_point_laplace,
                                                S = NULL,
                                                init = NULL,
                                                quiet = FALSE) {
  if (is.null(S)) {
    tic("Setting lower bound on row-wise residual variances")
    min_sd <- min(apply(dat, 1, sd))
    toc(quiet = quiet)
  }

  if (is.null(init)) {
    tic("Running PCA")
    init <- irlba::irlba(dat, nv = K)
    toc(quiet = quiet)
  }

  tic("Initializing grouped flash objects")
  group_idx <- split(seq_along(groups), groups)
  fl_list <- lapply(group_idx, function(idx) {
    grpdat <- dat[idx, ]
    init$u <- init$u[idx, ]
    fl <- flash_init(grpdat, var_type = 1, S = min_sd) |>
      flash_factors_init(init, ebnm_fn = ebnm_fn) |>
      flash_fit()
    fl$group_idx <- idx
    return(fl)
  })
  # Note that identities of groups can be retrieved via names(fl_list).
  toc(quiet = quiet)

  tic("Setting up cluster")
  cl <- makeCluster(nCores, type = "SOCK")
  zz <- clusterExport(
    cl,
    unclass(lsf.str(envir = asNamespace("flashier"), all = TRUE)),
    envir = as.environment(asNamespace("flashier"))
  )
  # Balance load based on number of features per group:
  grp_sizes <- sapply(fl_list, function(fl) nrow(fl$EF[[1]]))
  grp_order <- order(grp_sizes, decreasing = TRUE)
  assignments <- rep(list(NULL), nCores)
  for (grp in grp_order) {
      smallest_load <- which.min(
        sapply(assignments, function(core) sum(grp_sizes[core]))
      )
      assignments[[smallest_load]] <- c(assignments[[smallest_load]], grp)
  }
  for (core in seq_along(cl)) {
    sendCall(
      cl[[core]],
      assign,
      list("fl_list", fl_list[assignments[[core]]], envir = as.environment(1))
    )
  }
  zz <- checkForRemoteErrors(lapply(cl, recvResult))
  toc(quiet = quiet)

  if (!quiet) {
    cat("Cluster initialization complete. Data object may be removed from memory.")
  }

  return(cl)
}
