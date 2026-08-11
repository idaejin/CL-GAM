# Paths for CL-GAM experiments (Madrid MEDEA + local spclmm)
#
# Resolve relative to this file's directory so scripts work from any cwd.

.clgam_experiments_root <- local({
  # When sourced from R/, experiments root is parent; when from scripts/, grandparent
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f) || !nzchar(f)) {
    # fallback: assume working directory is experiments/ or experiments/scripts/
    wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    if (basename(wd) == "scripts" || basename(wd) == "R") {
      return(dirname(wd))
    }
    return(wd)
  }
  d <- dirname(normalizePath(f, winslash = "/", mustWork = TRUE))
  if (basename(d) %in% c("R", "scripts")) dirname(d) else d
})

CLGAM_ROOT <- normalizePath(
  file.path(.clgam_experiments_root, ".."),
  winslash = "/",
  mustWork = TRUE
)

CLGAM_EXP <- normalizePath(.clgam_experiments_root, winslash = "/", mustWork = TRUE)

CLGAM_SMIMR <- file.path(CLGAM_ROOT, "Diego Ayma", "SMiMR")
CLGAM_MADRID <- file.path(CLGAM_SMIMR, "Community of Madrid data analysis")
CLGAM_TIDY <- file.path(CLGAM_MADRID, "Tidy - Data")
CLGAM_CARTO <- file.path(CLGAM_MADRID, "Carto")
# Prefer experiment fork (no rgeos); fall back to archival SMiMR source
CLGAM_SPCLMM_SRC <- {
  fork <- file.path(CLGAM_EXP, "pkg", "spclmm")
  arch <- file.path(CLGAM_SMIMR, "spclmm")
  if (file.exists(file.path(fork, "DESCRIPTION"))) fork else arch
}
CLGAM_SPCLMM_ARCHIVAL <- file.path(CLGAM_SMIMR, "spclmm")
CLGAM_OUTPUT <- file.path(CLGAM_EXP, "output")
CLGAM_RLIBS <- file.path(CLGAM_EXP, "R_libs")

stopifnot(
  dir.exists(CLGAM_TIDY),
  dir.exists(CLGAM_CARTO),
  dir.exists(CLGAM_SPCLMM_SRC),
  file.exists(file.path(CLGAM_SPCLMM_SRC, "DESCRIPTION"))
)

if (!dir.exists(CLGAM_OUTPUT)) dir.create(CLGAM_OUTPUT, recursive = TRUE)

CLGAM_FAST <- identical(Sys.getenv("CLGAM_FAST", unset = "0"), "1")
CLGAM_NDX_SPATIAL <- if (CLGAM_FAST) c(8L, 8L) else c(20L, 20L)
CLGAM_NDX_NL <- if (CLGAM_FAST) 6L else 12L

message("CL-GAM experiments root: ", CLGAM_EXP)
message("Madrid tidy data:          ", CLGAM_TIDY)
message("spclmm source:             ", CLGAM_SPCLMM_SRC)
message("FAST mode:                 ", CLGAM_FAST, " (ndx spatial = ", paste(CLGAM_NDX_SPATIAL, collapse = ","), ")")
