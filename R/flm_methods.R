#' Print method for flm objects
#'
#' @description
#' Displays a concise summary of a fitted functional linear model object.
#'
#' @param x An object of class "flm", typically the result of a call to \code{\link{flm}}.
#' @param ... Additional arguments passed to model printing functions.
#'
#' @details
#' Prints basic information about the model including the call, number of observations,
#' response grid length, and a summary of each predictor.
#'
#' @return The function is called for printing model information.
#'
#' @seealso \code{\link{flm}}, \code{\link{summary.flm}}
#'
#' @examples
#' \dontrun{
#' # Fit functional factor regression model
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3))
#'
#' # Print basic model information
#' print(ffr_model)
#' }
#'
#' @rdname print.flm
#' @method print flm
#' @export
print.flm <- function(x, ...) {
  cat("\nFunction-on-Function Linear Regression\n\n")
  cat("Call:\n")
  print(x$info$call)
  cat("\nModel Information:\n")
  cat("Number of observations:", x$info$n_obs, "\n")
  cat("Response grid length:", x$info$grid_info$Y$length, "\n")
  cat("\nPredictors:\n")
  for (name in names(x$info$grid_info$X)) {
    pred_info <- x$info$grid_info$X[[name]]
    if (pred_info$type == "scalar") {
      cat(sprintf("%s: scalar predictor\n", name))
    } else {
      cat(sprintf(
        "%s: grid length = %d, factors = %d\n",
        name, pred_info$length, pred_info$K
      ))
    }
  }
}
#' Summary method for flm objects
#'
#' @description
#' Computes and returns comprehensive summary statistics for a fitted functional linear model.
#'
#' @param object An object of class "flm", typically the result of a call to \code{\link{flm}}.
#' @param ... Additional arguments passed to summary functions.
#'
#' @details
#' Calculates various goodness-of-fit measures including R-squared (both functional and average),
#' adjusted R-squared, RMSE, MAE, and model complexity metrics such as the total number of
#' factors and effective degrees of freedom.
#'
#' \subsection{Functional R-squared}{
#' For each point \eqn{r} in the response grid, the functional R-squared is calculated as:
#'
#' \deqn{R^2(r) = 1 - \frac{RSS(r)}{TSS(r)}}
#'
#' where \eqn{RSS(r)} is the residual sum of squares at grid point \eqn{r}:
#'
#' \deqn{RSS(r) = \sum_{t=1}^{T} (Y_t(r) - \widehat{Y}_t(r))^2}
#'
#' and \eqn{TSS(r)} is the total sum of squares at grid point \eqn{r}:
#'
#' \deqn{TSS(r) = \sum_{t=1}^{T} (Y_t(r) - \bar{Y}(r))^2}
#'
#' The average R-squared across all grid points is then computed as:
#'
#' \deqn{R^2 = \frac{1}{GL} \sum_{r=1}^{G} R^2(r)}
#'
#' where \eqn{G} is the number of points in the response grid.
#' }
#'
#' \subsection{Adjusted R-squared}{
#' The adjusted R-squared accounts for model complexity by penalizing for the number of factors
#' and scalar predictors:
#'
#' \deqn{R^2_{adj}(r) = 1 - \frac{RSS(r)/(T - p - q - 1)}{TSS(r)/(T - 1)}}
#'
#' where:
#' \itemize{
#'   \item \eqn{T} is the number of observations
#'   \item \eqn{p} is the total number of factors across all functional predictors
#'   \item \eqn{q} is the number of scalar predictors
#' }
#'
#' The effective degrees of freedom is calculated as \eqn{T - p - q - 1}.
#' }
#'
#' \subsection{Error Measures}{
#' The function also calculates:
#' \itemize{
#'   \item Root Mean Square Error (RMSE): \deqn{RMSE = \sqrt{\frac{1}{T \cdot G} \sum_{t=1}^{T} \sum_{r=1}^{G} (Y_t(r) - \widehat{Y}_t(r))^2}}
#'   \item Mean Absolute Error (MAE): \deqn{MAE = \frac{1}{T \cdot G} \sum_{t=1}^{T} \sum_{r=1}^{G} |Y_t(r) - \widehat{Y}_t(r)|}}
#' }
#'
#' @return An object of class "summary.flm" containing:
#'   \item{info}{Model information from the original flm object}
#'   \item{coefficients}{Regression coefficients from the original flm object}
#'   \item{inference}{Inference results if available}
#'   \item{residuals}{Model residuals}
#'   \item{stats}{List of model performance statistics including:
#'     \itemize{
#'       \item{R_squared}{Average R-squared across the response grid}
#'       \item{R_squared_func}{Functional R-squared values for each point in the response grid}
#'       \item{adj_R_squared}{Adjusted R-squared accounting for model complexity}
#'       \item{RMSE}{Root mean square error}
#'       \item{MAE}{Mean absolute error}
#'       \item{Y_mean_sq_error}{Mean squared error for the response}
#'       \item{X_mean_sq_errors}{Mean squared errors for each functional predictor's factor model}
#'       \item{total_factors}{Total number of factors used across all functional predictors}
#'       \item{n_scalar}{Number of scalar predictors}
#'       \item{effective_df}{Effective degrees of freedom}
#'     }
#'   }
#'
#' @seealso \code{\link{flm}}, \code{\link{print.summary.flm}}, \code{\link{plot.summary.flm}}
#'
#' @examples
#' \dontrun{
#' # Fit functional factor regression model
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3))
#'
#' # Compute summary statistics
#' model_summary <- summary(ffr_model)
#' }
#'
#' @rdname summary.flm
#' @method summary flm
#' @export
summary.flm <- function(object, ...) {
  # Calculate total sum of squares
  TSS <- colSums((object$fitted$Y + object$residuals$Y - pracma::repmat(object$info$grid_info$Y$mean, object$info$n_obs, 1))^2)
  # Calculate residual sum of squares
  RSS <- colSums(object$residuals$Y^2)
  # Calculate R-squared
  R_squared_func <- 1 - RSS / TSS
  R_squared_mean <- mean(R_squared_func)

  # Calculate total factors (only for functional predictors)
  total_factors <- sum(sapply(object$info$grid_info$X, function(x) {
    if (x$type == "functional") x$K else 0
  }))

  # Adjusted R-squared (adjusted for factors and scalar predictors)
  n_scalar <- sum(sapply(object$info$grid_info$X, function(x) x$type == "scalar"))
  adj_R_squared <- mean(1 - (RSS / (object$info$n_obs - total_factors - n_scalar)) /
    (TSS / (object$info$n_obs - 1)))

  # Calculate root mean square error
  RMSE <- sqrt(mean(object$residuals$Y^2))

  # Calculate mean absolute error
  MAE <- mean(abs(object$residuals$Y))

  # Create summary object
  sum_obj <- list(
    info = object$info,
    coefficients = object$coefficients,
    inference = object$inference,
    residuals = object$residuals,

    # Comprehensive model statistics
    stats = list(
      R_squared = R_squared_mean,
      R_squared_func = R_squared_func,
      adj_R_squared = adj_R_squared,
      RMSE = RMSE,
      MAE = MAE,
      Y_mean_sq_error = mean(object$residuals$Y^2),
      X_mean_sq_errors = lapply(object$residuals$X, function(x) mean(x^2)),
      total_factors = total_factors,
      n_scalar = n_scalar,
      effective_df = object$info$n_obs - total_factors - n_scalar - 1
    )
  )

  class(sum_obj) <- "summary.flm"
  return(sum_obj)
}
#' Print method for summary.flm objects
#'
#' @description
#' Displays a formatted summary of a functional linear model.
#'
#' @param x An object of class "summary.flm", typically the result of a call to \code{\link{summary.flm}}.
#' @param ... Additional arguments passed to summary printing functions.
#'
#' @details
#' Presents a detailed overview of model fit statistics, model complexity metrics,
#' and predictor information. Also provides guidance on how to visualize different
#' aspects of the model using the \code{\link{plot.summary.flm}} function.
#'
#' @return The function is called for printing summary information.
#'
#' @seealso \code{\link{flm}}, \code{\link{summary.flm}}, \code{\link{plot.summary.flm}}
#'
#' @examples
#' \dontrun{
#' # Fit functional factor regression model
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3))
#'
#' # Compute and print summary statistics
#' model_summary <- summary(ffr_model)
#' print(model_summary)
#' }
#'
#' @rdname print.summary.flm
#' @method print summary.flm
#' @export
#'
print.summary.flm <- function(x, ...) {
  cat("\nFunction-on-Function Linear Regression\n")
  cat("=======================================\n\n")
  cat("Call:\n")
  print(x$info$call)

  cat("\nModel Fit:\n")
  cat("------------\n")
  cat(sprintf("R-squared:          %7.4f\n", x$stats$R_squared))
  cat(sprintf("Adjusted R-squared: %7.4f\n", x$stats$adj_R_squared))
  cat(sprintf("RMSE:               %7.4f\n", x$stats$RMSE))
  cat(sprintf("MAE:                %7.4f\n", x$stats$MAE))

  cat("\nModel Complexity:\n")
  cat("----------------\n")
  cat(sprintf("Functional predictors factors: %d\n", x$stats$total_factors))
  cat(sprintf("Number of scalar predictors: %d\n", x$stats$n_scalar))
  cat(sprintf("Effective degrees of freedom: %d\n", x$stats$effective_df))

  cat("\nPredictor Information:\n")
  cat("--------------------\n")
  for (name in names(x$info$grid_info$X)) {
    pred_info <- x$info$grid_info$X[[name]]
    cat(sprintf("%s:\n", name))
    if (pred_info$type == "scalar") {
      cat("  Type: Scalar\n")
    } else {
      cat(sprintf("  Type: Functional\n"))
      cat(sprintf("  Grid length: %d\n", pred_info$length))
      cat(sprintf("  Factors (K): %d\n", pred_info$K))
    }
  }

  cat("\nVisualization Options:\n")
  cat("-------------------\n")
  cat("For functional predictors:\n")
  cat("  plot(summary(model), predictor = 'predictor_name', which = 'beta', conf.region = TRUE)  # Plot heatmap coefficient\n")
  cat("  plot(summary(model), predictor = 'predictor_name', which = 'beta_3D', conf.region = TRUE)  # Plot 3D coefficient\n")
  cat("  plot(summary(model), predictor = 'predictor_name', which = 't')     # Plot t-values\n")
  cat("  plot(summary(model), predictor = 'predictor_name', which = 'p')     # Plot p-values\n")
  cat("  plot(summary(model), which = 'R2')     # Plot functional R-squared\n")
  cat("\nFor scalar predictors:\n")
  cat("  plot(summary(model), predictor = 'intercept', which = 'beta')  # Plot intercept function\n")
  cat("  plot(summary(model), predictor = 'predictor_name', which = 'beta')  # Plot coefficient function\n")
}
