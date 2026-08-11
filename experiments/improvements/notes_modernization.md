# Modernisation notes (spclmm + paper scripts)

## Provenance
Factual: original PAPER Codes*.R use `maptools::readShapePoly`, absolute Windows `setwd`, and `pois_SAP` from local `spclmm` 0.1.0.

## Changes already in `experiments/`
1. **Paths** — `R/00_paths.R` points at SMiMR tidy data / Carto (no data copy).
2. **Cartography** — `sf::st_read` + `as(..., "Spatial")` for `choromap` (avoids retired `maptools` / `rgdal`).
3. **Package fork** — `pkg/spclmm` 0.1.8: `sparse.backend=` Matrix/spam/dense; partition `rowsum`; `sparse_P2D()`; Schur SOP. See `notes_bottleneck_backends.md`.
4. **FAST mode** — `CLGAM_FAST=1` shrinks `ndx` for dry runs (not paper numbers).
5. **Malone competitor** — Method 1 defaults (`poisson`+`round`+`spatial=none`). Bugfix: `t(C)` mass-balance. Stop on `crit_map` (report Codes3 `crit_fitgap` issue). Provenance: `improvements/MALONE_ADAPTATIONS.md`.

## Package improvements worth doing (NEXT)
Tracked here; do **not** silently rewrite scientific results.

| Issue | Where | Suggested fix |
|-------|--------|----------------|
| `class(Hinv) == "try-error"` | `pois_SAP` / related | use `inherits(Hinv, "try-error")` |
| Hard deps on spatial stack | DESCRIPTION / Imports | declare `sf`/`sp`/`classInt`/`RColorBrewer` explicitly; soft-deprecate maptools |
| Sparse algebra | composition `C` | accept `Matrix` and keep `C` sparse end-to-end |
| API docs | roxygen TODOs | document Case A/B/C arguments (`lcovfine`, `nlcovfine`, `C`) |
| Sex model | `pois_incat_SAP` | unit test against Codes1 expected AIC when available |
| Estimation | future | optional TMB Laplace backend (see `NOTES_ESTIMATION.md`) — not for V1 replication |

## Replication targets (ndx = 20,20)
From Codes2 comments:
- model0 AIC ≈ 374.58, ed ≈ 70.5 (~197 s historically)
- model1 AIC ≈ 396.92, ed ≈ 69.3 (~523 s)

FAST mode will **not** match these; use default `CLGAM_FAST=0` for numerical comparison.

## Do not
- Duplicate `Diego/` or zip trees into `experiments/`.
- Delete or move SMiMR data.
- Treat smoke-test AIC under FAST as manuscript evidence.
