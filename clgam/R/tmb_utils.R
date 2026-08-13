# Runtime TMB compile for method = "Laplace" (Suggests: TMB).

.clgam_tmb_env <- new.env(parent = emptyenv())

.clgam_tmb_src_dir <- function() {
  p <- system.file("tmb", "clgam_pois.cpp", package = "clgam")
  if (nzchar(p) && file.exists(p)) return(dirname(p))
  cand <- c(
    file.path(getwd(), "inst", "tmb"),
    file.path(getwd(), "clgam", "inst", "tmb")
  )
  for (d in cand) {
    if (file.exists(file.path(d, "clgam_pois.cpp"))) return(d)
  }
  stop("Cannot locate inst/tmb/clgam_pois.cpp.", call. = FALSE)
}

.clgam_tmb_compile_flags <- function() {
  cxx <- tryCatch(
    system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "config", "CXX17"),
      stdout = TRUE, stderr = FALSE
    ),
    error = function(e) character()
  )
  if (length(cxx) && any(grepl("clang|g\\+\\+|gcc", cxx, ignore.case = TRUE))) {
    return("-Wno-unused-but-set-variable")
  }
  ""
}

#' Compile and load the Laplace TMB template (once per session)
#' @keywords internal
.clgam_ensure_tmb_dll <- function(name = "clgam_pois", force = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop(
      "method='Laplace' requires package TMB. install.packages(\"TMB\")",
      call. = FALSE
    )
  }
  loaded <- vapply(getLoadedDLLs(), function(x) x[["name"]], character(1))
  if (!force && (name %in% loaded || isTRUE(.clgam_tmb_env[[name]]))) {
    return(invisible(name))
  }
  src_dir <- .clgam_tmb_src_dir()
  cpp <- file.path(src_dir, paste0(name, ".cpp"))
  if (!file.exists(cpp)) stop("TMB template not found: ", cpp, call. = FALSE)
  work <- file.path(tempdir(), "clgam_tmb")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  file.copy(cpp, file.path(work, basename(cpp)), overwrite = TRUE)
  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  setwd(work)
  TMB::compile(basename(cpp), flags = .clgam_tmb_compile_flags())
  dyn.load(TMB::dynlib(name))
  .clgam_tmb_env[[name]] <- TRUE
  invisible(name)
}
