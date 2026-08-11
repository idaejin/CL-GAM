#!/usr/bin/env Rscript
# Exploratory mgcv::gam at census-tract level (observed CT counts = oracle support).
# Not CL-GAMM: used to screen covariate shapes / collinearity before Case A.
#
# From experiments/:  Rscript scripts/04_mgcv_explore_covariates.R

root <- if (basename(getwd()) == "scripts") ".." else "."
source(file.path(root, "R/00_paths.R"))
source(file.path(root, "R/02_load_madrid.R"))

suppressPackageStartupMessages(library(mgcv))

dat <- clgam_load_madrid(load_maps = FALSE)
cov <- as.data.frame(dat$covariates)
for (nm in names(cov)) cov[[nm]] <- 100 * cov[[nm]]

df <- data.frame(
  y = dat$yc,
  e = dat$ec,
  lon = dat$xxc[, 1],
  lat = dat$xxc[, 2],
  cov
)
df <- df[df$e > 0, ]
message("n = ", nrow(df))

vars <- c("manual", "unemployed", "salaried", "pollution", "ageing", "low_edu")

m0 <- gam(
  y ~ offset(log(e)) + s(lon, lat, k = 100),
  family = poisson(), data = df, method = "REML"
)

uni <- lapply(setNames(vars, vars), function(v) {
  message("uni: ", v)
  fml <- as.formula(sprintf(
    "y ~ offset(log(e)) + s(lon, lat, k = 80) + s(%s, k = 10)", v
  ))
  gam(fml, family = poisson(), data = df, method = "REML")
})

m_paper <- gam(
  y ~ offset(log(e)) + s(lon, lat, k = 80) +
    s(unemployed, k = 10) + s(ageing, k = 10),
  family = poisson(), data = df, method = "REML"
)

m_lean <- gam(
  y ~ offset(log(e)) + s(lon, lat, k = 80) +
    s(unemployed, k = 10) + s(ageing, k = 10) + s(pollution, k = 10),
  family = poisson(), data = df, method = "REML"
)

m_ses <- gam(
  y ~ offset(log(e)) + s(lon, lat, k = 80) +
    manual + s(unemployed, k = 10) + s(ageing, k = 10) +
    s(pollution, k = 10) + s(low_edu, k = 10),
  family = poisson(), data = df, method = "REML"
)

summarize_fit <- function(m, name) {
  sm <- summary(m)
  list(
    name = name,
    aic = AIC(m),
    bic = BIC(m),
    deviance_expl = sm$dev.expl,
    r_sq = sm$r.sq,
    edf_total = sum(sm$edf),
    s_table = as.data.frame(sm$s.table)
  )
}

fits <- c(
  list(summarize_fit(m0, "spatial_only")),
  lapply(vars, function(v) summarize_fit(uni[[v]], paste0("spatial+", v))),
  list(
    summarize_fit(m_paper, "spatial+unemp+ageing"),
    summarize_fit(m_lean, "spatial+unemp+ageing+pollution"),
    summarize_fit(m_ses, "spatial+SES_block")
  )
)

cmp <- do.call(rbind, lapply(fits, function(f) {
  data.frame(
    model = f$name,
    aic = f$aic,
    bic = f$bic,
    dev_expl = f$deviance_expl,
    r_sq = f$r_sq,
    edf = f$edf_total,
    stringsAsFactors = FALSE
  )
}))
cmp <- cmp[order(cmp$aic), ]
cmp$delta_aic <- cmp$aic - min(cmp$aic)

partial_curve <- function(m, vn, n = 80) {
  newdata <- data.frame(
    e = 1,
    lon = mean(df$lon),
    lat = mean(df$lat),
    manual = mean(df$manual),
    unemployed = mean(df$unemployed),
    salaried = mean(df$salaried),
    pollution = mean(df$pollution),
    ageing = mean(df$ageing),
    low_edu = mean(df$low_edu)
  )
  newdata <- newdata[rep(1, n), ]
  rownames(newdata) <- NULL
  xr <- range(df[[vn]], na.rm = TRUE)
  newdata[[vn]] <- seq(xr[1], xr[2], length.out = n)
  pr <- predict(m, newdata = newdata, type = "terms", se.fit = TRUE)
  hit <- grep(paste0("^s\\(", vn, "\\)"), colnames(pr$fit), value = TRUE)
  data.frame(
    x = newdata[[vn]],
    fit = as.numeric(pr$fit[, hit]),
    se = as.numeric(pr$se.fit[, hit]),
    variable = vn
  )
}

curves_lean <- rbind(
  transform(partial_curve(m_lean, "unemployed"), model = "spatial+unemp+ageing+pollution"),
  transform(partial_curve(m_lean, "ageing"), model = "spatial+unemp+ageing+pollution"),
  transform(partial_curve(m_lean, "pollution"), model = "spatial+unemp+ageing+pollution")
)

out <- list(
  note = paste(
    "Exploratory mgcv::gam at CT with observed counts (oracle).",
    "Poisson + offset(log(e)); covariates as %."
  ),
  n = nrow(df),
  correlation = cor(df[, vars]),
  comparison = cmp,
  curves_lean = curves_lean,
  s_table_lean = summary(m_lean)$s.table,
  s_table_paper = summary(m_paper)$s.table,
  s_table_ses = summary(m_ses)$s.table
)

saveRDS(out, file.path(CLGAM_OUTPUT, "mgcv_explore.rds"))
write.csv(cmp, file.path(CLGAM_OUTPUT, "mgcv_compare.csv"), row.names = FALSE)
write.csv(curves_lean, file.path(CLGAM_OUTPUT, "mgcv_curves_lean.csv"), row.names = FALSE)
write.csv(out$correlation, file.path(CLGAM_OUTPUT, "mgcv_cor.csv"))

print(cmp)
message("Wrote ", file.path(CLGAM_OUTPUT, "mgcv_explore.rds"))
