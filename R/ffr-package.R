#' Functional Factor Regression
#'
#' @description
#' A comprehensive framework for estimation, prediction and statistical inference in linear models with functional response
#' and functional/scalar predictors. The package implements the functional factor regression approach for function-on-function
#' linear regression models.
#'
#'
#' @details
#' The \pkg{ffr} package provides tools for regression modeling of functional data through a factor-based
#' approach that enables statistical inference. Key features include:
#'
#' \itemize{
#'   \item \strong{Model Estimation:} Fit function-on-function linear regression models with the \code{\link{flm}} function
#'   \item \strong{Factor Determination:} Consistently estimate the number of factors using the eigenvalue
#'         difference approach with \code{\link{fed}} function
#'   \item \strong{Parameter Tuning:} Optimize model hyperparameter with cross-validation using \code{\link{tune.fed}}
#'   \item \strong{Prediction:} Generate forecasts on new data with \code{\link{predict.flm}}
#'   \item \strong{Visualization:} Examine model components through specialized plotting methods,
#'         including coefficient surfaces, t-values, and confidence regions
#'   \item \strong{Data Preprocessing:} Transform sparsely or irregularly sampled functional data to a consistent
#'         grid with \code{\link{fpreprocess}}
#' }
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{flm}}: Fit a function-on-function linear regression model
#'   \item \code{\link{fed}}: Estimate the number of factors for functional predictors
#'   \item \code{\link{tune.fed}}: Tune the gamma parameter for factor estimation
#'   \item \code{\link{fpreprocess}}: Preprocess functional data to a consistent grid
#'   \item \code{\link{predict.flm}}: Generate predictions from fitted models
#'   \item \code{\link{summary.flm}}: Summarize model fit statistics
#'   \item \code{\link{print.summary.flm}}: Comprehensively print model summary
#'   \item \code{\link{plot.summary.flm}}: Create visualizations of model components
#' }
#'
#' @references
#' Otto, S., & Winter, L. (in press). Functional factor regression with an application to electricity price curve modeling.
#'
#' Wu, J. (2018). Eigenvalue difference test for the number of common factors in the approximate factor models. Economics Letters, 169, 63–67.
#'
#' @keywords internal
"_PACKAGE"
#' @name ffr
#' @import stats
#' @importFrom pracma repmat
#' @importFrom graphics abline contour image legend matplot par
#' @importFrom grDevices colorRampPalette
NULL
