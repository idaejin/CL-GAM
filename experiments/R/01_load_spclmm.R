# Ensure local spclmm is available (source under Diego Ayma/SMiMR/spclmm).

if (!exists("CLGAM_SPCLMM_SRC")) {
  paths <- c("R/00_paths.R", "../R/00_paths.R", "00_paths.R")
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stop("Source R/00_paths.R first (from experiments/).")
  source(hit)
}

.clgam_local_lib <- function() {
  lib <- file.path(CLGAM_EXP, "R_libs")
  if (!dir.exists(lib)) dir.create(lib, recursive = TRUE)
  .libPaths(c(lib, .libPaths()))
  lib
}

.clgam_ensure_spclmm <- function(force_reinstall = FALSE) {
  lib <- .clgam_local_lib()
  need <- isTRUE(force_reinstall) || !requireNamespace("spclmm", quietly = TRUE)
  desc <- read.dcf(file.path(CLGAM_SPCLMM_SRC, "DESCRIPTION"))
  local_ver <- unname(desc[1, "Version"])
  # Rebuild if choromap/source is newer than installed package
  if (!need) {
    inst_dir <- file.path(lib, "spclmm")
    src_mtime <- max(file.mtime(c(
      list.files(file.path(CLGAM_SPCLMM_SRC, "R"), full.names = TRUE, pattern = "[.]R$"),
      list.files(file.path(CLGAM_SPCLMM_SRC, "src"), full.names = TRUE, pattern = "[.](cpp|h|c)$")
    )), na.rm = TRUE)
    inst_mtime <- if (dir.exists(inst_dir)) {
      file.mtime(file.path(inst_dir, "DESCRIPTION"))
    } else {
      as.POSIXct(0)
    }
    if (is.finite(src_mtime) && src_mtime > inst_mtime) {
      message("spclmm source newer than install; reinstalling…")
      need <- TRUE
    }
  }
  if (!need) {
    inst_ver <- as.character(utils::packageVersion("spclmm"))
    if (!identical(local_ver, inst_ver)) {
      message(
        "Installed spclmm ", inst_ver, " != source ", local_ver,
        "; reinstalling from ", CLGAM_SPCLMM_SRC
      )
      need <- TRUE
    }
  }
  if (need) {
    message("Installing spclmm into ", lib, " …")
    utils::install.packages(
      CLGAM_SPCLMM_SRC,
      lib = lib,
      repos = NULL,
      type = "source"
    )
  }
  suppressPackageStartupMessages(library(spclmm, lib.loc = lib))
  invisible(TRUE)
}

.clgam_ensure_spclmm()
