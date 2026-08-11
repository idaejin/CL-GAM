# CL-GAMM estimation notes (SOP canonical)

**Status:** DECISION for revive V1 — estimation engine is **SOP / SAP** (Rodríguez-Álvarez et al., 2015), not numerical REML maximisation and not TMB.

**Code truth:** [`Diego Ayma/SMiMR/spclmm/R/pois_SAP.R`](Diego%20Ayma/SMiMR/spclmm/R/pois_SAP.R) (Case A); sex contrast via `pois_incat_SAP.R`. Parallel `pois_REML.R` is legacy / cross-check only.

**Manuscript:** [`manuscript-CL-GAMM-V1/CL-GAMM-V1.tex`](manuscript-CL-GAMM-V1/CL-GAMM-V1.tex)

---

## Model (Case A)

\[
\boldsymbol{\mu}=\mathbf{C}\bigl(\boldsymbol{e}_f\odot\exp(\mathbf{X}_f\boldsymbol{\beta}+\mathbf{Z}_f\boldsymbol{\alpha})\bigr),
\quad
\boldsymbol{\alpha}\sim\mathcal{N}(\boldsymbol{0},\mathbf{G}),
\quad
\mathbf{G}^{-1}=\sum_{d=1}^{2+K}\tau_d^{-2}\boldsymbol{\Lambda}_d.
\]

- \(d=1,2\): anisotropic spatial P-spline penalties (longitude / latitude).
- \(d=3,\ldots,2+K\): one variance component per nonlinear fine-scale smooth \(g_k\).
- Linear fine covariates enter \(\mathbf{X}_f\) only (no \(\tau^2\)).

In `pois_SAP`, stored `la` / `var.comp` are the **\(\tau^2\)** values; penalties use `Ginv` entries proportional to \(1/\tau_d^2\).

---

## Algorithm 1 — PIRLS + SOP (matches `pois_SAP`)

**Input:** \(\boldsymbol{y}\) (coarse), \(\mathbf{C}\), \(\boldsymbol{e}_f\), fine coordinates / covariates, bases \(\mathbf{X}_f,\mathbf{Z}_f\), penalty structure \(\{\boldsymbol{\Lambda}_d\}\).

1. Initialise \(\boldsymbol{\beta},\boldsymbol{\alpha}\) (e.g. zero) and \(\boldsymbol{\tau}^{2(0)}\) (e.g. ones). Set \(\boldsymbol{\eta}=\mathbf{X}_f\boldsymbol{\beta}+\mathbf{Z}_f\boldsymbol{\alpha}\), \(\boldsymbol{\gamma}=\boldsymbol{e}_f\odot\exp(\boldsymbol{\eta})\), \(\boldsymbol{\mu}=\mathbf{C}\boldsymbol{\gamma}\).

2. **Outer PIRLS loop** until \(\|\boldsymbol{\eta}^{\mathrm{new}}-\boldsymbol{\eta}\|^2/\|\boldsymbol{\eta}\|^2<\varepsilon_{\eta}\):

   a. Working vector (composite-link Poisson linearisation)
   \[
   \boldsymbol{z}
   =
   \bigl(\mathrm{diag}(\boldsymbol{\mu})^{-1}\mathbf{C}\,\mathrm{diag}(\boldsymbol{\gamma})\bigr)\boldsymbol{\eta}
   +
   \mathrm{diag}(\boldsymbol{\mu})^{-1}(\boldsymbol{y}-\boldsymbol{\mu}).
   \]

   b. Working cross-products via `clmm_mat` (GLAM-friendly):  
      \(\breve{\mathbf{X}}'\mathbf{W}\breve{\mathbf{X}}\), \(\breve{\mathbf{X}}'\mathbf{W}\breve{\mathbf{Z}}\), \(\breve{\mathbf{Z}}'\mathbf{W}\breve{\mathbf{Z}}\), \(\breve{\mathbf{X}}'\mathbf{W}\boldsymbol{z}\), \(\breve{\mathbf{Z}}'\mathbf{W}\boldsymbol{z}\) with \(\mathbf{W}=\mathrm{diag}(\boldsymbol{\mu})\).

   c. **Inner SOP loop** until mean \(|\tau_d^{2\,\mathrm{new}}-\tau_d^{2}|<\varepsilon_{\tau}\):

      - Form \(\mathbf{G}^{-1}=\sum_d \tau_d^{-2}\boldsymbol{\Lambda}_d\) (block structure as in Lee–Durbán / Ayma et al.).
      - Solve the mixed system for \(\boldsymbol{\beta},\boldsymbol{\alpha}\) (PQL / Henderson).
      - For each \(d\), update (Schall / SAP)
        \[
        \widehat{\mathrm{ED}}_d
        =
        \mathrm{tr}\!\bigl(
          \text{diag contribution of }\boldsymbol{\Lambda}_d\text{ to }\mathbf{Z}'\mathbf{PZ}
        \bigr),
        \qquad
        \tau_d^{2\,\mathrm{new}}
        =
        \frac{\boldsymbol{\alpha}'\boldsymbol{\Lambda}_d\boldsymbol{\alpha}}{\widehat{\mathrm{ED}}_d}.
        \]
        (Implemented in `pois_SAP` via diagonal extracts `dZtNZ` and masks `G1inv.n`, `G2inv.n`, `Gkinv.n`.)

   d. Update \(\boldsymbol{\eta},\boldsymbol{\gamma},\boldsymbol{\mu}\) from new \(\boldsymbol{\beta},\boldsymbol{\alpha}\).

3. **Output:** \(\widehat{\boldsymbol{\eta}}\) at fine scale, \(\widehat{\boldsymbol{\tau}}^2\), EDs, optional AIC/BIC from total ED, pointwise SE for \(\boldsymbol{\eta}\) and \(g_k\).

---

## Equivalence REML ↔ SOP

For fixed working weights, SOP updates are the fixed-point equations of REML for variance components in the linear mixed model on \(\boldsymbol{z}\) (Schall 1991; Rodríguez-Álvarez et al. 2015). Outer PIRLS updates the working model for the composite-link Poisson mean. Numerical maximisation of the REML profile (`pois_REML`) targets the same fixed point but is slower / less stable with many \(\tau^2\).

---

## Madrid API consistency *(Factual)*

| Script | Estimator | Role |
|---|---|---|
| `PAPER - Codes2.R` | **`pois_SAP`** only | ATA maps + Case A covariates |
| `PAPER - Codes1.R` | **`pois_incat_SAP`** | Sex contrast |
| `PAPER - Simulations*.R` | **`pois_SAP`** | Competitors vs Malone / multistep |
| — | `pois_REML` / `pois_incat_REML` | **Not used** in paper scripts |

Canonical claim in text = what Madrid ran.

---

## NEXT (not V1)

- **Camarda & Durbán** ([arXiv:2412.04956](https://arxiv.org/abs/2412.04956)): latent working response + Kronecker \(P\) — implemented as `pois_PCLM_latent` in experiment fork (`experiments/improvements/notes_camarda_2412.04956.md`). Cited in `manuscript-CL-GAMM-V1` as complementary computation (array/separable \(C\)); **not** the V1 estimation claim (still SOP).
- **Experiment fork speed *(Factual)*:** `experiments/pkg/spclmm` **0.2.0** compiles Schur / aggregation / Gram kernels via RcppArmadillo; Madrid mun / mun→ct ≈ **8–10×** vs pure R with identical paper AIC (`scripts/15_rcpp_bench.R`). Prefer this over sparse-basis reformulation for irregular areal \(C\).
- **Evaluation *(DECISION)*:** sims → primary **MSE/MAE(\(\hat\eta\) vs \(\eta_{\mathrm{true}}\))** + optional cover95 vs truth; real ATA/models → coarse deviance/log-score, mass \(\max\|C\hat\mu-y\|\), AIC only within the same engine; fine oracle secondary. Details: `experiments/improvements/notes_metrics_ATA.md` · Brain `methods/ATA-evaluation-metrics.md`. Helper: `experiments/R/08_metrics.R`.
- **Public package *(DECISION)*:** new tree [`clgam/`](clgam/) (`Package: clgam`); simulated ATA via `simulate_ata()`; no MEDEA in-repo. See `experiments/improvements/notes_clgam_package.md`.
- **TMB + Laplace** for the same GLMM (marginal likelihood / SEs) as robustness appendix or software track with [[TMB-P-splines]].
- Do not redefine the paper’s novelty as “TMB composite link” in V1.
