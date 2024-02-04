#' @importFrom tictoc tic toc
#' @importFrom snow clusterEvalQ
#' @export
flash_impute_grouped_data <- function(cl,
                                      quiet = FALSE) {
  tic("Imputing missing data")
  zz <- clusterEvalQ(cl, {
    fl_list <<- mapply(function(fl, which_na) {
      new_data <- get.Y(fl)
      new_data[which_na] <- fitted(fl)[which_na]
      return(update_data_workhorse(fl, new_data, wrapup = FALSE))
    }, fl_list, na_list, SIMPLIFY = FALSE)
    NULL
  })
  toc(quiet = quiet)

  return(invisible(NULL))
}
