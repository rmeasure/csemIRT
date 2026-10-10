# csemIRT

**Bootstrap-validated conditional standard errors of measurement for polytomous IRT models**

`csemIRT` estimates the conditional standard error of measurement (CSEM) and conditional reliability for the graded response model (GRM), generalized partial credit model (GPCM), partial credit model (PCM), graded rating scale model (GRSM) and rating scale model (RSM). It extends `csemGRM` 0.1.0, which supported the GRM only.

Main features:

- **Five polytomous models on one metric.** All CSEM profiles are reported on a standardized latent metric (unit variance). This makes models that estimate the latent variance (PCM, RSM) directly comparable with models that fix it (GRM, GPCM, GRSM).
- **Uncertainty of the CSEM profile.** You can choose a parametric or nonparametric bootstrap band (`csem_bootstrap()`) or a fast delta-method band (`csem_delta()`). Failed replications are counted and reported by reason.
- **Model comparison.** `csem_compare()` refits all models to the same bootstrap resamples. It reports, pointwise, whether the choice of model changes the estimated precision.
- **Precision regions.** `csem_regions()` translates conditional-reliability targets into CSEM thresholds. It reports the theta intervals and population shares that meet each target, with a conservative band-based version.
- **Plots.** `csem_plot()` draws ggplot2 graphics for every result.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rmeasure/csemIRT")
```

## Quick start

```r
library(csemIRT)

items <- attitude_data[, 1:10]

fit  <- csem_estimate(items, model = "graded")    # or "gpcm", "pcm", "grsm", "rsm"
csem_analytical(fit)                              # CSEM = 1 / sqrt(I(theta))
csem_reliability(fit)                             # conditional reliability
csem_regions(fit)                                 # where is rho >= .90/.80/.70?

dlt  <- csem_delta(fit)                           # delta-method band (seconds)
boot <- csem_bootstrap(fit, n_boot = 1000,        # bootstrap band
                       type = "parametric", cores = 2)
csem_plot(boot)

# Does model choice change the precision conclusions?
cmp <- csem_compare(wellbeing_data[, 1:12],
                    models = c("graded", "gpcm", "pcm"),
                    n_boot = 1000, cores = 2)
cmp
csem_plot(cmp, type = "band")
csem_plot(cmp, type = "difference")
```

## Functions

| Function | Purpose |
|---|---|
| `csem_estimate()` | Fit GRM / GPCM / PCM / GRSM / RSM (via `mirt`); information, person estimates, marginal reliability, AIC/BIC |
| `csem_analytical()` | Analytical CSEM profile |
| `csem_delta()` | Delta-method confidence band (log-scale Wald) |
| `csem_bootstrap()` | Parametric or nonparametric bootstrap band, with failure diagnostics |
| `csem_reliability()` | Conditional reliability (`"ratio"` or `"complement"` definition), with bands |
| `csem_compare()` | Multi-model comparison with a paired bootstrap of CSEM differences |
| `csem_regions()` | Theta regions and population coverage meeting reliability targets |
| `csem_plot()` | Plots: `csem`, `info`, `items`, `reliability`, `band`, `difference` |

## Built-in data

- `attitude_data`: 10 items, 5 categories, N = 500, generated from the GRM.
- `character_data`: 12 items, 5 categories, N = 400, generated from the GRM.
- `wellbeing_data`: 12 items, 4 categories, N = 600, generated from the GPCM with heterogeneous slopes and a ceiling effect.

## Dependencies

R (>= 4.1), `mirt` (>= 1.40), `ggplot2` (>= 3.4).

## License

MIT
