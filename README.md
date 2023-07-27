# flashier

Functions documentation and vignettes are available [here](https://willwerscheid.github.io/flashier/).

Please note that function names changed in Version 0.2.44. If you are using `flashier` code from previous versions, please make the following substitutions, and note that many parameter names have changed as well:

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
