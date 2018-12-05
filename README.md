# flashier

This package is in active development. The interface is likely to change from time to time.

Please note that `flashier` currently uses a development branch of package `ebnm`. To get the default `point.normal` priors (and some other prior types) working, you will need to run:

```devtools::install_github("stephenslab/ebnm, ref = "regular-normal")```

When using the more flexible `ashr` priors (such as `normal.mixture` and `nonnegative`), I recommend that you include `ebnm.ash = list(optmethod = "mixSQP")` in your `flashier` calls if you care about runtime. To use this option, you will need to install `mixsqp`. I strongly recommend using the current version:

```devtools::install_github("stephenslab/mixsqp")```
