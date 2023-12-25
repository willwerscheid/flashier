flash_recover_fl_object_from_cluster <- function(cl) {
  fl_list <- clusterEvalQ(cl, fl_list)
  fl_list <- Reduce(c, fl_list)

  group_idx <- lapply(fl_list, `[[`, "group_idx")
  idx <- order(Reduce(c, group_idx))

  groups <- character(length(idx))
  for (i in seq_along(group_idx)) {
    groups[group_idx[[i]]] <- names(group_idx)[i]
  }

  fl <- fl_list[[1]]
  fl$Y <- Reduce(rbind, lapply(fl_list, `[[`, "Y"))[idx, ]
  if (fl$Z != 1) {
    fl$Z <- Reduce(rbind, lapply(fl_list, `[[`, "Z"))[idx, ]
  }
  fl$n.nonmissing <- Reduce(c, lapply(fl_list, `[[`, "n.nonmissing"))[idx]
  fl$Y2 <- Reduce(c, lapply(fl_list, `[[`, "Y2"))[idx]
  fl$R2 <- Reduce(c, lapply(fl_list, `[[`, "R2"))[idx]
  fl$est.tau <- Reduce(c, lapply(fl_list, `[[`, "est.tau"))[idx]
  fl$tau <- Reduce(c, lapply(fl_list, `[[`, "est.tau"))[idx]
  fl$obj <- clusterEvalQ(cl, elbo)[[1]]
  fl$EF[[1]] <- Reduce(rbind, lapply(fl_list, function(grp_fl) grp_fl$EF[[1]]))[idx, ]
  fl$EF2[[1]] <- Reduce(rbind, lapply(fl_list, function(grp_fl) grp_fl$EF2[[1]]))[idx, ]
  fl$KL[[1]] <- rowSums(sapply(fl_list, function(grp_fl) grp_fl$KL[[1]]))
  fl$KL[[2]] <- clusterEvalQ(cl, KL_F)[[1]]
  g_L <- lapply(fl_list, function(fl) lapply(fl$g, `[[`, 1))
  for (k in seq_along(fl$g)) {
      fl$g[[k]][[1]] <- lapply(g_L, `[[`, k)
      prior_fam <- ebnm:::infer_prior_family(fl_list[[1]]$g[[k]][[1]])
      fl$ebnm.fn[[k]][[1]] <- flash_ebnm(group = groups,prior_family = prior_fam)
  }
  fl$is.zero <- apply(sapply(fl_list, function(grp_fl) grp_fl$is.zero), 1, all)
  fl$is.valid <- apply(sapply(fl_list, function(grp_fl) grp_fl$is.valid), 1, all)

  return(fl)
}
