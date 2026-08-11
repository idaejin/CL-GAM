# Malone adaptations (provenance)

**Status:** DECISION for competitor code in `experiments/R/03_malone.R`.  
**Claim labels:** Factual = in Malone 2012 or CL-GAMM Method 1 text/Codes3; Adaptation = our/disease-mapping departure from Malone 2012; Bugfix = Codes3 implementation error only.

## From Malone et al. (2012) — keep

- Iterative loop: fine-scale regression on covariates → multiplicative rescale so fine values **aggregate to the coarse map** → repeat until stable.
- Use of fine-resolution covariates to drive disaggregation.

## CL-GAMM Method 1 / Codes3 — disease-mapping adaptations *(not in Malone 2012)*

These are **already in the CL-GAMM draft** (Method 1), not invented in the 2026 revive:

| Adaptation | Why |
|---|---|
| Irregular nested areal \(\mathbf{C}\) (mun→CT) | Malone 2012 is raster coarse→fine |
| Init \(\hat{\mathbf{y}}_f = \mathbf{C}^{+}\mathbf{y}\) (equal-split counts) | Malone pastes the coarse *value* by nearest neighbour (continuous); counts cannot be pasted or totals explode |
| Poisson log-link + `offset(log e_f)` | Malone target is continuous earth-resource attribute |
| `round()` working counts before GAM | Needed for Poisson integer response in Method 1 / Codes3 |
| Stopping tolerance \(0.001\) on mean absolute change | As written in Method 1 |

## Bugfix only *(not a scientific novelty)*

- Codes3: `yc <- yc * rep(w, rowSums(C))` — correct **only if** CT are contiguous by municipality (true for Madrid `C.m`, fragile in general).
- Revive: `yc * as.vector(t(C) %*% w)` — same math when ordered; safe if reordered.

## Optional flags in `clgam_malone_fit` — report if used

| Flag | Default Method 1 | If changed, report as |
|---|---|---|
| `spatial = "none"` | yes | `"additive"` = Codes3 Madrid extra; `"te"` = further revive option |
| `family = "poisson"` | yes | `"quasipoisson"` + no round = revive numerical adaptation |
| `round_counts = TRUE` | yes | `FALSE` only with quasipoisson |
| `init = "equal"` | yes | `"exposure"` = mun SMR × \(e_f\) (revive intermediate) |
| `deliver = "fitted"` | yes | `"balanced"` = post-rescale map with exact mass (revive intermediate) |
| `stop_on = "map"` | yes | `"balanced_map"` / `"mass"` = revive stopping rules |

## Revive intermediates (`scripts/07_malone_intermediates.R`) *(not Malone 2012, not Method 1)*

Heuristic-only tweaks — **still no \(\mathbf{C}\) in the likelihood**:

1. **Deliver balanced map** — after the last GAM, report \(\mathbf{y}^{\mathrm{bal}}\) with \(\mathbf{C}\mathbf{y}^{\mathrm{bal}}=\mathbf{y}\) (pattern from GAM, totals forced).
2. **Quasipoisson + no round** — continuous working response (avoids integer noise).
3. **Exposure init** — \(y_j=e_j\cdot(y_i/e_i)\) instead of equal-split counts.
4. **Stop on balanced-map / mass** — alternative convergence monitors.

Note: multiplicative rescale of \(\boldsymbol{\mu}\) is algebraically the same as rate rescale when \(\mu_j=e_j r_j\); the distinct exposure idea is in **init** (and deliver), not a second rescale formula.

## Stopping criterion *(report)*

- **Malone / Method 1 text & Codes3:** stop when \((1/m)\sum|\hat y_j^{l}-\hat y_j^{l-1}|\le 0.001\). Codes3 coded this as `mean(|y_bal - fitted|)` after each rescale+GAM.
- **Observed with Poisson+round:** that gap (`crit_fitgap`) stays large because the GAM is not mass-constrained; Codes3 likely hit `max_iter` more often than true 0.001 convergence.
- **Revive default stop:** `mean(|μ_t - μ_{t-1}|)` on successive GAM fitted maps (`crit_map`). All three diagnostics are stored. This is an **implementation clarification**, not a change to the Malone rescale step.

## Do not claim

- That Method 1 **is** Malone 2012 without the disease adaptations above.
- That `quasipoisson` / `s(lon,lat)` are from Malone 2012.
- That mass-balance via `t(C)` changes the Madrid numerical answer when CT stay in mun order (bugfix / robustness only).
- That Codes3’s printed `crit < 0.001` on `|y_bal - fitted|` was typically met under Poisson+round.
