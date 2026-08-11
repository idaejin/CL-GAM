#' Composite Link Generalized Additive (Mixed) Models
#'
#' Fits CL-GAM / CL-GAMM models for nested areal count disaggregation with
#' P-splines and separation of overlapping precision matrices (SOP).
#'
#' @keywords internal
#' @useDynLib clgam, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom stats AIC BIC coef cor deviance dpois fitted logLik nobs
#'   predict residuals runif rpois sd median
#' @importFrom graphics plot points abline lines par rect text mtext polygon
#'   legend
#' @importFrom grDevices hcl.colors adjustcolor dev.interactive devAskNewPage
#' @import Matrix
"_PACKAGE"
