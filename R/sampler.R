build.sampler <- function(flash) {
  # Beware of unfulfilled promise leak.
  force(flash)

  return(function(nsamp) {
    # Get samples as list of dimensions with sublists of factors.
    samp <- rapply(all.post.samplers(flash),
                   function(f) do.call(f, list(nsamp = nsamp)),
                   how = "list")
    # Re-organize the list so that each element corresponds to a single sample.
    return(lapply(1:nsamp, function(trial) {
      lapply(1:get.dim(flash),
             function(n) do.call(cbind,
                                 lapply(samp[[n]], function(k) k[trial, ])))
    }))
  })
}

all.post.samplers <- function(flash) {
  return(lapply(1:get.dim(flash),
                function(n) lapply(1:get.n.factors(flash),
                                   function(k) one.post.sampler(flash, k, n))))
}

one.post.sampler <- function(flash, k, n) {
  factor <- extract.factor(flash, k)
  if (all.fixed(factor, n)) {
    sampler <- function(nsamp) {matrix(get.EF(factor)[[n]],
                                       nrow = nsamp,
                                       ncol = get.dims(flash)[n],
                                       byrow = TRUE)}
  } else if (is.zero(factor)) {
    sampler <- function(nsamp) {matrix(0,
                                       nrow = nsamp,
                                       ncol = get.dims(flash)[n])}
  }
  else {
    ebnm.res <- solve.ebnm(factor, n, flash, output = "posterior_sampler")
    sampler <- ebnm.res$posterior_sampler
  }
  return(sampler)
}
