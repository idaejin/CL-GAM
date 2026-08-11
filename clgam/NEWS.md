# clgam 0.1.9

* Rename the PIRLS+SOP engine to `pois_SOP` / `pois_incat_SOP`
  (separation of overlapping precision matrices). Legacy names
  `pois_SAP` / `pois_incat_SAP` remain as aliases.
* Rename Rcpp Schur kernel to `sop_solve_schur_cpp`.
* Document SOP estimation in DESCRIPTION and README.
* Optional `pois_TMB()` Laplace check (Suggests: TMB).
* Generate Rd help pages via roxygen2.

# clgam 0.1.8

* Previous release (SOP kernels under historical SAP names).
