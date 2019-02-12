F.sampler <- function(flash) {
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
                                   function(k) post.sampler(flash, k, n))))
}

post.sampler <- function(flash, k, n) {
  factor <- extract.factor(flash, k)
  if (is.zero(factor)) {
    sampler <- function(nsamp) {matrix(0,
                                       nrow = nsamp,
                                       ncol = get.dims(flash)[n])}
  } else if (all.fixed(factor, n)) {
    sampler <- function(nsamp) {matrix(get.fix.vals(flash, k),
                                       nrow = nsamp,
                                       ncol = get.dims(flash)[n],
                                       byrow = TRUE)}
  } else {
    ebnm.res <- solve.ebnm(factor, n, flash, return.sampler = TRUE)
    sampler <- ebnm.res$post_sampler
  }
  return(sampler)
}
