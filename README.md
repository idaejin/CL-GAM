# CL-GAM (R package `clgam`)

This repository contains:

- The public R package implementation in `./clgam/`.
- Replication / experiments utilities in `./experiments/` (Madrid, pennLC, sims, competitors).

For standard package usage (API, evaluation hierarchy, install details), see:
`./clgam/README.md`.

## Install (source)

Install the package from the repo (recommended):

```r
# from the repository root
install.packages("clgam", repos = NULL, type = "source")
```

If you already `cd` into `clgam/`, then:

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick check (sanity)

```r
library(clgam)
dat <- simulate_ata(n_coarse = 8, n_fine_per = 4, seed = 1, include_covariate = TRUE)
fit <- clgam(
  y = dat$y, coords = cbind(dat$x1, dat$x2), C = dat$C, exposure = dat$efine,
  smooth = dat$nlcovfine, knots = c(6, 6), knots_nl = 8
)
plot(fit, which = 4, g_true = dat$g_true)
```

## Notes for paper-ready evaluation

Simulation and model comparison metrics are documented in:
`experiments/improvements/notes_metrics_ATA.md` (Research Brain: `ATA-evaluation-metrics`).

## License

GPL (>= 2)

