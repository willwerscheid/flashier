# TODO: timeout parameter?
#' @importFrom tictoc tic toc
#' @importFrom snow clusterEvalQ clusterExport
#' @export
#'
flash_backfit_grouped_data <- function(cl,
                                       maxiter,
                                       tol = 1e-6,
                                       tolparam = c("factors", "elbo"),
                                       extrapolate = TRUE,
                                       warmstart = TRUE,
                                       quiet = FALSE) {
  tolparam <- match.arg(tolparam)

  tic("Backfitting flash object")
  Fdim <- clusterEvalQ(cl, dim(fl_list[[1]]$EF[[2]]))[[1]]
  kset <- seq_len(Fdim[2])
  ebnm_fn <- clusterEvalQ(cl, fl_list[[1]]$ebnm.fn)[[1]][[1]][[2]]

  clusterExport(cl, "warmstart", envir = environment())
  zz <- clusterEvalQ(cl, {
    fl_list <<- lapply(fl_list, set.warmstart, warmstart)
    NULL
  })

  iter <- 0
  is_converged <- FALSE

  if (extrapolate) {
    extrapolate.control <- getOption("extrapolate.control", list())
    extrapolate.param <- set.extrapolate.param(extrapolate.control)
    extrapolate.param <- init.beta(extrapolate.param)
    clusterExport(cl, "extrapolate.param", envir = environment())
    zz <- clusterEvalQ(cl, {old_fl_list <<- fl_list; NULL})
  }

  elbo <- -Inf
  new_F <- matrix(0, nrow = Fdim[1], ncol = Fdim[2])

  while (iter < maxiter && !is_converged) {
    iter <- iter + 1

    old_elbo <- elbo
    old_F <- new_F

    if (extrapolate) {
      zz <- clusterEvalQ(cl, {
        pre_list <<- fl_list
        fl_list <<- mapply(
          extrapolate.f, fl_list, old_fl_list,
          SIMPLIFY = FALSE, MoreArgs = list(par = extrapolate.param)
        )
        NULL
      })

      update_res <- update_factors_in_grouped_data(cl, kset, new_F, ebnm_fn, extrapolated = TRUE)
      is_extrapolation_successful <- update_res$elbo > old_elbo

      if (is_extrapolation_successful) {
        extrapolate.param <- accelerate(extrapolate.param)
      } else {
        zz <- clusterEvalQ(cl, {fl_list <<- pre_list; NULL})
        extrapolate.param <- decelerate(extrapolate.param)
      }
      clusterExport(cl, "extrapolate.param", envir = environment())
    }

    if (!extrapolate || !is_extrapolation_successful) {
      update_res <- update_factors_in_grouped_data(cl, kset, new_F, ebnm_fn, extrapolated = FALSE)
    }

    new_F <- update_res$new_F
    KL_F <- update_res$KL_F
    elbo <- update_res$elbo

    elbo_diff <- update_res$elbo - old_elbo
    F_chg <- calc.max.abs.chg(new_F, old_F)

    if (tolparam == "elbo") {
      is_converged <- elbo_diff < tol
    } else { # tolparam == "factors"
      is_converged <- F_chg < tol
    }


    if (!quiet) {
      cat(paste0(
        formatC(iter, width = 4), ".",
        " Fchg = ", formatC(F_chg, format = "e", digits = 3),
        ", ELBO = ", formatC(elbo, format = "f", digits = 2), "\n"
      ))
    }
  }
  toc(quiet = quiet)

  clusterExport(cl, c("KL_F", "elbo"), envir = environment())
  return(invisible(new_F))
}

update_factors_in_grouped_data <- function(cl, kset, new_F, ebnm_fn, extrapolated) {
  KL_F <- rep(0, length(kset))

  for (k in kset) {
    clusterExport(cl, "k", envir = environment())

    args_list <- clusterEvalQ(cl, {
      # Update L.
      fct_list <<- lapply(fl_list, function(fl) {
        fct <- extract.factor(fl, k)
        if (!is.zero(fct)) {
          return(update.factor.one.n(fct, 1, fl))
        } else {
          return(fct)
        }
      })
      # Return arguments to EBNM problem for F.
      mapply(function(fct, fl) {
        if (is.zero(fct)) {
          args <- NULL
        } else {
          args <- calc.ebnm.args(fct, n = 2, fl, include.fixed = FALSE)
          if (any(args$s == 0)) {
            args <- NULL
          }
        }
        return(args)
      }, fct_list, fl_list, SIMPLIFY = FALSE)
    })

    # Combine arguments across groups and solve EBNM problem for F.
    args_list <- Reduce(c, args_list)
    args_list <- args_list[!sapply(args_list, is.null)]
    if (length(args_list) > 0) {
      tauEL2 <- 1 / sapply(args_list, `[[`, "s")^2
      s <- 1 / sqrt(rowSums(tauEL2))
      x <- rowSums(sapply(args_list, `[[`, "x") * tauEL2) * s^2
      ebnm_res <- ebnm_fn(
        x, s, output = c(
          "posterior_mean", "posterior_second_moment", "fitted_g", "log_likelihood"
        )
      )
      KL_F[k] <- ebnm_res$log_likelihood - normal.means.loglik(
        x, s, ebnm_res$posterior$mean, ebnm_res$posterior$second_moment
      )
    } else {
      ebnm_res <- list(
        posterior = data.frame(mean = 0, second_moment = 0),
        fitted_g = NULL
      )
      KL_F[k] <- 0
    }

    # Update F and return ELBO.
    clusterExport(cl, c("ebnm_res", "extrapolated"), envir = environment())
    elbo_list <- clusterEvalQ(cl, {
      fct_list <- lapply(fct_list, function(fct) {
        if (!is.null(fct)) {
          fct <- set.EF(fct, ebnm_res$posterior$mean, 2)
          fct <- set.EF2(fct, ebnm_res$posterior$second_moment, 2)
          fct <- set.g(fct, ebnm_res$fitted_g, 2)
        }
        return(fct)
      })

      if (extrapolated) {
        fct_list <- lapply(fct_list, set.to.valid)
      } else {
        fct_list <- mapply(
          update.R2.tau.and.obj, fct_list, fl_list,
          SIMPLIFY = FALSE
        )
      }

      fl_list <<- mapply(
        insert.factor, fl_list, fct_list,
        SIMPLIFY = FALSE, MoreArgs = list(update.tau = !extrapolated)
      )

      lapply(fl_list, get.obj)
    })

    # Keep track of F.
    new_F[, k] <- ebnm_res$posterior$mean
  }

  if (extrapolated) {
    elbo_list <- clusterEvalQ(cl, {
      fl_list <- lapply(fl_list, init.tau)
      fl_list <<- lapply(fl_list, function(fl) set.obj(fl, calc.obj(fl)))
      lapply(fl_list, get.obj)
    })
  }

  elbo <- sum(unlist(Reduce(c, elbo_list))) + sum(KL_F)

  return(list(new_F = new_F, KL_F = KL_F, elbo = elbo))
}
