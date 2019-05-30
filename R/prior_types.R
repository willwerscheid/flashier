# Prior types

normal <- function(...) {
  param <- list(prior_type = "normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

point.normal <- function(...) {
  param <- list(prior_type = "point_normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

point.laplace <- function(...) {
  param <- list(prior_type = "point_laplace")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

nonzero.mode <- function(...) {
  param <- list(prior_type = "point_normal", fix_mu = FALSE)
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.pn,
                   ebnm.param = modifyList(param, list(...)))))
}

normal.mix <- function(...) {
  param <- list(mixcompdist = "normal")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

uniform.mix <- function(...) {
  param <- list(mixcompdist = "uniform")
  return(list(list(sign = 0,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

nonnegative <- function(...) {
  param <- list(mixcompdist = "+uniform")
  return(list(list(sign = 1,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}

nonpositive <- function(...) {
  param <- list(mixcompdist = "-uniform")
  return(list(list(sign = -1,
                   ebnm.fn = ebnm.ash,
                   ebnm.param = modifyList(param, list(...)))))
}
