# flashier

Documentation and vignettes are available [here][pkgdown-site].

## Quick start

Install the latest version of flashier from GitHub:

```r
install.packages("remotes")
remotes::install_github("willwerscheid/flashier")
```

Once you have installed the package, load the package in R:

```r
library(flashier)
```

Next, run an example analysis of the GTEx data set:

```r
data(gtex)
fl <- flash(gtex,greedy_Kmax = 3,backfit = TRUE)
```

For a more detailed introduction to flashier, see the
[introductory vignette][pkgdown-vignette-intro] and, later, the
[other vignettes][pkgdown-vignettes].

To learn more, visit the [package website][pkgdown-site], and view the
"flash" help page in R:

```r
help("flash")
```

## Citing this work

If you find the flashier package or any of the source code in this
repository useful for your work, please cite:

> Jason Willwerscheid, Peter Carbonetto and Matthew Stephens. ebnm: an
> R package for solving the empirical Bayes normal means problem using
> a variety of prior families. [arXiv:2110.00152][ebnm-preprint].

## License

Copyright (c) 2018-2024, Jason Willwerscheid, Peter Carbonetto and
Matthew Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].

## Notes

Please note that function names changed in flashier version 0.2.44. If
you are using flashier code from previous versions, please make the
following substitutions, and note that many parameter names have
changed as well:

Old Name	| New Name
--- | ---
as.ebnm.fn | flash_ebnm
conv.crit.elbo	| flash_conv_crit_elbo_diff
conv.crit.factors	| flash_conv_crit_max_chg_F
conv.crit.loadings	| flash_conv_crit_max_chg_L
display.elbo	| flash_verbose_elbo
display.elbo.diff	| flash_verbose_elbo_diff
display.F.max.chg	| flash_verbose_max_chg_F
display.L.max.chg	| flash_verbose_max_chg_L
display.max.chg	| flash_verbose_max_chg
ff.elbo	| flash_fit_get_elbo
ff.est.tau	| flash_fit_get_est_tau
ff.fixed.tau	| flash_fit_get_fixed_tau
ff.g	| flash_fit_get_g
ff.KL	| flash_fit_get_KL
ff.p2m	| flash_fit_get_p2m
ff.pm	| flash_fit_get_pm
ff.tau	| flash_fit_get_tau
flash.add.greedy	| flash_greedy
flash.backfit	| flash_backfit
flash.fit	| flash_fit
flash.fix.factors	| flash_factors_fix
flash.init	| flash_init
flash.init.factors	| flash_factors_init
flash.nullcheck	| flash_nullcheck
flash.remove.factors	| flash_factors_remove
flash.reorder.factors	| flash_factors_reorder
flash.set.factors.to.zero	| flash_factors_set_to_zero
flash.set.verbose	| flash_set_verbose
init.fn.default	| flash_greedy_init_default
init.fn.irlba	| flash_greedy_init_irlba
init.fn.softImpute	| flash_greedy_init_softImpute

## Credits

The flashier R package was developed by Jason Willwerscheid, Peter
Carbonetto and Matthew Stephens, with many other contributors.

[mit-license]: https://opensource.org/licenses/mit-license.html
[ebnm-preprint]: https://arxiv.org/abs/2110.00152
[pkgdown-site]: https://willwerscheid.github.io/flashier/
[pkgdown-vignettes]: https://willwerscheid.github.io/flashier/articles/index.html
[pkgdown-vignette-intro]: https://willwerscheid.github.io/flashier/articles/flashier_intro.html
