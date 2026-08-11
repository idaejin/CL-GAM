# Runtime TMB compile/load for clgam (DLL separate from Rcpp useDynLib).

.tmb_env <- new.env(parent = emptyenv())

.tmb_src_dir <- function() {
  cand <- system.file("TMB", package = "clgam")
  if (nzchar(cand) && file.exists(file.path(cand, "clgam_pois.cpp"))) {
    return(cand)
  }
  pkg <- tryCatch(find.package("clgam", quiet = TRUE), error = function(e) "")
  if (nzchar(pkg)) {
    for (d in c(
      file.path(pkg, "inst", "TMB"),
      file.path(pkg, "TMB"),
      file.path(dirname(pkg), "inst", "TMB")
    )) {
      if (file.exists(file.path(d, "clgam_pois.cpp"))) return(d)
    }
  }
  for (d in c(
    file.path(getwd(), "inst", "TMB"),
    file.path(getwd(), "clgam", "inst", "TMB"),
    file.path(getwd(), "TMB")
  )) {
    if (file.exists(file.path(d, "clgam_pois.cpp"))) return(d)
  }
  stop("Cannot locate TMB templates for clgam (inst/TMB/clgam_pois.cpp).", call. = FALSE)
}

.tmb_compile_flags <- function() {
  cxx <- tryCatch(
    system2(file.path(R.home("bin"), "R"),
            c("CMD", "config", "CXX17"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(cxx) && any(grepl("clang|g\\+\\+|gcc", cxx, ignore.case = TRUE))) {
    return("-Wno-unused-but-set-variable")
  }
  ""
}

#' Ensure a TMB DLL is compiled and loaded
#' @param name Template basename without extension (e.g. \code{"clgam_pois"})
#' @param force Recompile even if already loaded
#' @keywords internal
ensure_tmb_dll <- function(name = "clgam_pois", force = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required for pois_TMB. install.packages(\"TMB\")",
         call. = FALSE)
  }
  loaded <- vapply(getLoadedDLLs(), function(x) x[["name"]], character(1))
  if (!force && (name %in% loaded || isTRUE(.tmb_env[[name]]))) {
    return(invisible(name))
  }

  src_dir <- .tmb_src_dir()
  cpp <- file.path(src_dir, paste0(name, ".cpp"))
  if (!file.exists(cpp)) stop("TMB template not found: ", cpp, call. = FALSE)

  work <- file.path(tempdir(), "clgam_tmb")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)
  file.copy(cpp, file.path(work, basename(cpp)), overwrite = TRUE)

  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  setwd(work)

  TMB::compile(basename(cpp), flags = .tmb_compile_flags())
  dyn.load(TMB::dynlib(name))
  .tmb_env[[name]] <- TRUE
  invisible(name)
}
