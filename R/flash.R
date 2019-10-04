flash <- function(data = NULL,
                  S = NULL,
                  prior.family = prior.point.normal(),
                  var.type = 0L,
                  backfit = FALSE,
                  greedy.Kmax = 50L,
                  verbose = 1L) {
  fl <- flash.init(data, var.type, S)

  fl <- flash.add.greedy(fl, greedy.Kmax, prior.family,
                         verbose.lvl = verbose,
                         output.lvl = 3L * !backfit)

  if (backfit) {
    fl <- flash.backfit(fl, verbose.lvl = verbose)
  }

  return(fl)
}
