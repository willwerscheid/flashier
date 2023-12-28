flash_impute_grouped_data <- function(cl,
                                      quiet = FALSE) {
  tic("Imputing flash object")
  zz <- clusterEvalQ(cl, {
    fl_list <<- mapply(function(fl, which_na) {
      new_data <- get.Y(fl)
      new_data[which_na] <- fitted(fl)[which_na]
      fl |> flash_update_data(new_data) |> flash_fit()
    }, fl_list, na_list, SIMPLIFY = FALSE)
    NULL
  })
  toc(quiet = quiet)

  return(NULL)
}
