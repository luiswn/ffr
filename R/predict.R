#' Predict method for flm objects
#'
#' @description
#' Generates predictions from a fitted functional linear model for new data.
#'
#' @param object An object of class "flm", typically the result of a call to \code{\link{flm}}.
#' @param newdata A list containing the predictors for which to generate predictions. Must contain
#'        all predictors used in the original model with the same format and grid dimensions.
#'
#' @details
#' This function applies a fitted functional linear model to new data to generate functional predictions.
#' For each functional predictor, \code{predict} calculates factor scores for the new observations.
#' The final prediction is computed by multiplying the assembled design matrix (including
#' the intercept, scalar predictors, and factor scores) by the estimated
#' coefficient matrix from the original model.
#'
#' @note
#' In order to compute the factor scores, the functional regressors will be centered using the means
#' from the training data, not the means of the new data, to ensure consistency.
#'
#' @return A matrix of predicted values with dimensions:
#'   \itemize{
#'     \item Rows: Number of new observations in \code{newdata}
#'     \item Columns: Number of grid points in the functional response
#'   }
#'   If the original response had column names, these are preserved in the predictions.
#'
#' @seealso \code{\link{flm}}, \code{\link{summary.flm}}
#'
#' @examples
#' \dontrun{
#' # Fit functional factor regression model
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3))
#'
#' # Generate new data
#' # Here simplified for technical illustrative purposes only. Consult the flm
#' # documentation example for a relevant workable DGP.
#' new_data <- list(
#'   X = matrix(rnorm(100 * 200), nrow = 100, ncol = 200), # Functional predictor
#'   w = rnorm(100) # Scalar predictor
#' )
#'
#' # Generate predictions
#' predictions <- predict(ffr_model, newdata = new_data)
#' }
#'
#' @rdname predict.flm
#' @method predict flm
#' @export
predict.flm <- function(object, newdata) {
  # Validate input
  if (!inherits(object, "flm")) {
    stop("'object' must be a flm model")
  }
  if (!is.list(newdata)) {
    stop("'newdata' must be a list containing observations")
  }

  # Extract predictor names from model
  predictor_vars <- names(object$info$grid_info$X)

  # Check if all required predictors are present in newdata
  missing_vars <- setdiff(predictor_vars, names(newdata))
  if (length(missing_vars) > 0) {
    stop("Missing predictors in newdata: ", paste(missing_vars, collapse = ", "))
  }

  # Separate scalar and functional predictors
  scalar_vars <- names(which(sapply(object$info$grid_info$X, function(x) x$type == "scalar")))
  functional_vars <- setdiff(predictor_vars, scalar_vars)

  # Get number of new observations from functional predictors
  if (length(functional_vars) > 0) {
    n_new <- unique(sapply(newdata[functional_vars], nrow))
  } else {
    # If only scalar predictors, get length of first predictor
    n_new <- length(newdata[[predictor_vars[1]]])
  }

  if (length(n_new) != 1) {
    stop("All predictors must have the same number of observations")
  }

  # Function to compute new factors for a functional predictor
  compute_new_factors <- function(X_new, var) {
    # Get grid length and mean from training
    X_gridlength <- object$info$grid_info$X[[var]]$length
    X_mean <- object$info$grid_info$X[[var]]$mean

    # Validate dimensions of new data
    if (ncol(X_new) != X_gridlength) {
      stop(sprintf(
        "Predictor %s has incorrect grid length. Expected %d, got %d",
        var, X_gridlength, ncol(X_new)
      ))
    }

    # Center new data using training mean
    X_centered <- sweep(X_new, 2, X_mean, "-")

    # Get eigenvectors from training
    psi <- object$fitted$factors$Psi[[var]]

    # Compute new factors
    f_new <- X_centered %*% psi / X_gridlength

    return(f_new)
  }

  # Create new design matrix z_hat starting with intercept
  new_z <- matrix(1, nrow = n_new, ncol = 1) # intercept column

  # Add scalar predictors directly
  for (var in scalar_vars) {
    new_z <- cbind(new_z, as.matrix(newdata[[var]]))
  }

  # Compute and add factors for functional predictors
  if (length(functional_vars) > 0) {
    new_factors <- lapply(functional_vars, function(var) {
      compute_new_factors(newdata[[var]], var)
    })

    for (f in new_factors) {
      new_z <- cbind(new_z, f)
    }
  }

  # Make predictions using B coefficients
  predictions <- new_z %*% object$coefficients$B

  # Add column names if they exist in the original model
  if (!is.null(colnames(object$fitted$Y))) {
    colnames(predictions) <- colnames(object$fitted$Y)
  }

  return(predictions)
}
