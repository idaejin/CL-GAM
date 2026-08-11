#!/usr/bin/env Rscript
# Sex contrast (PAPER Codes1): pois_incat_SAP, mun→ct, female vs male.
# Paper runtime ~15 min with ndx=20,20 (here typically a few minutes).
#
# From experiments/:  Rscript scripts/03_replicate_sex_contrast.R
#
# Note: pois_incat_SAP returns var.comp (not tau2); aic/bic are length-2 (per sex).

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/01_load_spclmm.R"))
source(file.path(root, "R/02_load_madrid.R"))

dat <- clgam_load_madrid(load_maps = TRUE)

n_ct <- ncol(dat$C_m)
ym_all <- c(dat$mun_female$observed, dat$mun_male$observed)
ec_all <- c(dat$ct_female$expected, dat$ct_male$expected)
sex_ata <- factor(c(rep("female", n_ct), rep("male", n_ct)))
coords <- rbind(dat$xxc, dat$xxc)

message("Female mun deaths: ", sum(dat$mun_female$observed),
        " | Male mun deaths: ", sum(dat$mun_male$observed))
message("Stacked y length (2 * n_mun): ", length(ym_all))
message("Fine stacked length: ", length(ec_all), " | ndx = ", paste(CLGAM_NDX_SPATIAL, collapse = ","))

message("=== pois_incat_SAP sex contrast (Codes1) ===")
fit <- pois_incat_SAP(
  y = ym_all,
  x1 = coords[, 1],
  x2 = coords[, 2],
  efine = ec_all,
  cat = sex_ata,
  Ccat1 = dat$C_m,
  Ccat2 = dat$C_m,
  ndx = CLGAM_NDX_SPATIAL,
  elements = TRUE,
  trace = TRUE
)

vc <- as.numeric(fit$var.comp)
message(sprintf(
  "AIC (F,M)=%s  BIC (F,M)=%s  ed=%.3f  edf=%s",
  paste(round(fit$aic, 3), collapse = ", "),
  paste(round(fit$bic, 3), collapse = ", "),
  fit$ed,
  paste(round(fit$edf, 3), collapse = ", ")
))
message("var.comp: ", paste(round(vc, 6), collapse = ", "))
message("Paper var.comp ref: 0.2718913, 0.08473693, 0.06023534, 0.1038672")

dif <- fit$eta[1:n_ct] - fit$eta[(n_ct + 1):(2 * n_ct)]
message(sprintf("dif range: [%.4f, %.4f] (paper ~ [-0.787, 0.422])",
                min(dif), max(dif)))

brks <- as.numeric(quantile(dif, seq(0, 1, 0.1)))
brks[1] <- brks[1] - 0.001
brks[length(brks)] <- brks[length(brks)] + 0.001

pdf(file.path(CLGAM_OUTPUT, "sex_difference_surface.pdf"))
choromap(
  sf = dat$map_ct$sf, values = dif, breaks = brks,
  border = "grey30", lwd = 0.1, title = "Female - Male eta (CT)"
)
dev.off()

if (!is.null(fit$sd.dif)) {
  sd_dif <- fit$sd.dif
  brks_sd <- seq(min(sd_dif), max(sd_dif), length.out = 11)
  brks_sd[1] <- brks_sd[1] - 0.001
  brks_sd[length(brks_sd)] <- brks_sd[length(brks_sd)] + 0.001
  pdf(file.path(CLGAM_OUTPUT, "sex_sd_difference_surface.pdf"))
  choromap(
    sf = dat$map_ct$sf, values = sd_dif, breaks = brks_sd,
    name = "Reds", n = 9, border = "grey30", lwd = 0.1,
    title = "SD of difference surface"
  )
  dev.off()
}

saveRDS(
  list(
    ndx = CLGAM_NDX_SPATIAL,
    fast = CLGAM_FAST,
    aic = fit$aic,
    bic = fit$bic,
    ed = fit$ed,
    edf = fit$edf,
    var.comp = vc,
    eta = fit$eta,
    dif = dif,
    sd.dif = fit$sd.dif,
    paper_var.comp = c(0.2718913, 0.08473693, 0.06023534, 0.1038672),
    paper_dif_range = c(-0.7870461, 0.4219904)
  ),
  file.path(CLGAM_OUTPUT, "sex_contrast_fit.rds")
)

message("Wrote sex contrast outputs under ", CLGAM_OUTPUT)
