# flashier

This package is in active development. The interface is likely to change from time to time.

Please note that `flashier` currently uses a development branch of package `ebnm`. To get the default `point.normal` priors (and some other prior types) working, you will need to run:

```devtools::install_github("stephenslab/ebnm", ref = "regular-normal")```

To get more flexible priors such as `normal.mixture` and `nonnegative` working, you will need to install `ashr`:

```install.packages("ashr")```

When using the more flexible `ashr` priors, I recommend that you include `ebnm.ash = list(optmethod = "mixSQP")` in your `flashier` calls if you care about runtime. To use this option, you will need to install `mixsqp`. I strongly recommend using the current version:

```devtools::install_github("stephenslab/mixsqp")```

Finally, I recommend taking a look at the vignettes before diving in. Be forewarned, however, that they take a few minutes to build (it took me a little less than 5 minutes on a 2015 MacBook Pro). If you choose to build them, make sure you've installed the above packages first. Then run:

```devtools::install_github("willwerscheid/flashier", build_vignettes = TRUE)```
