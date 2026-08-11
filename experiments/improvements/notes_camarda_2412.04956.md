# Camarda & Durbán (arXiv:2412.04956) — ideas for CL-GAM

**Paper:** Camarda, C. G. & Durbán, M. (2025). *Fast Estimation of the Composite Link Model for Multidimensional Grouped Counts.* [arXiv:2412.04956](https://arxiv.org/abs/2412.04956).

**Provenance:** Factual (paper); PRELIMINARY (Madrid port).

## Ideas adopted

1. **Latent working response** \(\tilde y = \gamma \odot (C^\top (y \oslash \mu))\) — redistribute coarse counts to fine support without building \(\breve B = W^{-1} C \Gamma B\).
2. **Kronecker P-spline penalty** \(P=\lambda_1(I\otimes D_1'D_1)+\lambda_2(D_2'D_2\otimes I)\) via existing `sparse_P2D()`.
3. **Preserve Poisson likelihood on aggregated \(y\)** (not EM on latent Poisson).

## What we did *not* port yet

- Full **GLAM** \(\rho/\mathcal{G}\) array algebra for rectangular \(C=C_2\otimes C_1\) (mortality age×year). Madrid ATA uses **irregular** mun→CT `C`.
- Automatic \(\lambda\) selection / SAP equivalence.
- Analytic SEs from \(V=(B^\top\Gamma C^\top W^{-1}C\Gamma B+P)^{-1}\) (paper §2.3).

## Code

- `pkg/spclmm/R/pois_PCLM_latent.R` — `pois_PCLM_latent()`
- `scripts/14_camarda_latent_pclm.R` — Madrid mun→ct vs `pois_SAP`

## Relation to SOP revive

**DECISION (unchanged):** manuscript V1 estimation claim remains **SOP/`pois_SAP`**.  
**NEXT:** Camarda latent engine as optional fast backend / competitor; λ calibration and GLAM path for separable ST arrays.
