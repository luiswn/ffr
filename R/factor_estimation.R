#' Eigenvalue difference estimator to determine the number of factors in functional factor regressions
#'
#' @description
#' Consistently estimates the correct number of factors for each functional predictor in a function-on-function
#' regression model using an eigenvalue difference (ED) approach.
#'
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model. Uses the same
#'        format as in \code{\link{flm}}. A formula specification \code{response ~ .} implies
#'        using all variables in \code{data} except the response.
#' @param data A list containing the variables in the model. Each named object within the list is considered
#'        either a functional or a scalar variable, following the same format as in \code{\link{flm}}.
#' @param gamma A positive numeric tuning parameter that controls the sensitivity of factor detection.
#'        Default is 1.
#' @param plot Logical indicating whether to plot the eigenvalue difference function for visual inspection.
#'        Default is FALSE.
#'
#' @details
#' The ED estimator for functional factor regressions as introduced by Otto & Winter (2025) determines the number of factors by
#' finding the largest distance between two subsequent transformed eigenvalues of the integral operator \deqn{D}, which is the
#' product of the cross-covariance operator between the regressor and the regressand with its adjoint. Consult the original paper
#' for details about the method.
#'
#' @return A list with components:
#'   \item{gamma}{The tuning parameter used}
#'   \item{K}{A named vector containing the estimated number of factors for each functional predictor}
#'
#' @seealso \code{\link{flm}} for fitting functional linear models, \code{\link{tune.fed}} for tuning the gamma parameter
#'
#' @examples
#' \dontrun{
#' # Simulate functional data
#'
#' # Define a Fourier basis function
#' fourier.basis <- function(evalgrid, nbasis) {
#'   basis <- matrix(NA, length(evalgrid), nbasis)
#'   for (k in 1:nbasis) {
#'     if (k == 1) {
#'       basis[, k] <- rep(1, length(evalgrid))
#'     } else {
#'       if (k %% 2 == 0) {
#'         basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * evalgrid)
#'       } else {
#'         basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * evalgrid)
#'       }
#'     }
#'   }
#'   return(basis)
#' }
#'
#' # Function to generate random bivariate beta coefficients
#' rand.biv.beta <- function(Y.grid, X.grid, K, smoothness_penalty = 1, seed) {
#'   Y.basis <- matrix(NA, length(Y.grid), K)
#'   X.basis <- matrix(NA, length(X.grid), K)
#'   for (k in 1:K) {
#'     if (k == 1) {
#'       Y.basis[, k] <- rep(1, length(Y.grid))
#'       X.basis[, k] <- rep(1, length(X.grid))
#'     } else {
#'       if (k %% 2 == 0) {
#'         Y.basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * Y.grid)
#'         X.basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * X.grid)
#'       } else {
#'         Y.basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * Y.grid)
#'         X.basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * X.grid)
#'       }
#'     }
#'   }
#'
#'   set.seed(seed)
#'
#'   weights <- exp(-seq(0, 2, length.out = K))
#'
#'   # Create coefficient matrix with rank K
#'   A <- matrix(rnorm(K * K), K, K)
#'
#'   # Ensure coefficient matrix has exactly rank K through SVD
#'   svd_A <- svd(A)
#'   A_rank_K <- svd_A$u %*% diag(weights) %*% t(svd_A$v)
#'
#'   # Apply smoothness penalty
#'   A_smooth <- A_rank_K * smoothness_penalty
#'
#'   biv.func <- Y.basis %*% A_smooth %*% t(X.basis)
#'
#'   return(biv.func)
#' }
#'
#' # Set up simulation parameters
#' Y.gridlength <- 200
#' X.gridlength <- 200
#' Y.range <- c(0, ((Y.gridlength - 1) / Y.gridlength))
#' X.range <- c(0, ((X.gridlength - 1) / X.gridlength))
#' Y.grid <- seq(Y.range[1], Y.range[2], length.out = Y.gridlength)
#' X.grid <- seq(X.range[1], X.range[2], length.out = X.gridlength)
#'
#' K <- 3
#' T <- 100 # Number of observations
#'
#' # Generate data
#' four.basis <- fourier.basis(X.grid, 3 * K)
#'
#' f <- matrix(rnorm(n = T * K, 0, 1), T, K)
#' z <- cbind(rep(1, nrow(f)), f)
#' psi <- four.basis[, 1:K]
#' pred.func <- f %*% t(psi)
#' eps <- matrix(rnorm(n = T * 2 * K, 0, 1), T, 2 * K)
#' eps.func <- eps %*% t(four.basis[, (K + 1):(3 * K)])
#'
#' # Create functional predictor
#' X <- pred.func + eps.func
#'
#' # Generate coefficient surface
#' beta <- rand.biv.beta(Y.grid, X.grid, K, seed = 1)
#'
#' intercept <- rep(3, Y.gridlength)
#'
#' # Generate error term
#' u <- matrix(rnorm(n = T * 2 * K, 0, 1), T, 2 * K)
#' u.func <- u %*% t(fourier.basis(Y.grid, 2 * K))
#'
#' # Generate response
#' Y <- intercept + (X %*% beta / X.gridlength) + u.func
#'
#'
#' # Estimate number of factors
#' data <- list(Y = Y, X = X)
#' ED_result <- fed(Y ~ X, data, gamma = 1)
#' print(ED_result$K) # correct would be K=3
#' }
#'
#' @references
#' Otto, S., & Winter, L. (2025). Functional factor regression with an application to electricity price curve modeling.
#' Wu, J. (2018). Eigenvalue difference test for the number of common factors in the approximate factor models. Economics Letters, 169, 63–67.
#'
#' @rdname fed
#' @importFrom pracma repmat
#' @importFrom utils tail
#' @export
fed <- function(formula, data, gamma = 1, plot = FALSE) {
  # Validate basic inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.list(data)) {
    stop("'data' must be a list containing observations")
  }

  # Extract variable names from formula
  vars <- all.vars(formula)
  response_var <- all.vars(formula[[2]]) # Left-hand side

  # Handle dot notation in formula
  if (length(all.vars(formula[[3]])) == 1 && all.vars(formula[[3]]) == ".") {
    # If "." is used, include all variables except the response
    predictor_vars <- setdiff(names(data), response_var)
  } else {
    # Otherwise use the variables specified in the formula
    predictor_vars <- setdiff(vars, response_var)
  }

  # Validate that all variables exist in data
  missing_vars <- setdiff(predictor_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  is_scalar <- function(var_data) {
    if (is.list(var_data)) {
      # Check if it's a data list object with explicit type
      if (!is.null(var_data$type)) {
        return(var_data$type == "scalar")
      }
    }
    # Check if it's a data frame with a single column
    if (ncol(as.data.frame(var_data)) == 1) {
      return(TRUE)
    }
    # Default to functional if no type specified and not a single-column list
    return(FALSE)
  }

  # Identify scalar and functional predictors
  scalar_vars <- predictor_vars[sapply(data[predictor_vars], is_scalar)]
  functional_vars <- setdiff(predictor_vars, scalar_vars)

  # Extract response and predictor data
  Y <- if (is_scalar(data[[response_var]])) data[[response_var]]$data else data[[response_var]]

  # Extract and validate functional predictors
  X_functional <- lapply(functional_vars, function(var) {
    if (is_scalar(data[[var]])) data[[var]]$data else data[[var]]
  })
  names(X_functional) <- functional_vars

  # Validate dimensions
  T <- if (is.matrix(Y)) nrow(Y) else length(Y)

  # Validate functional predictor dimensions
  if (!all(sapply(X_functional, nrow) == T)) {
    stop("All functional variables must have the same number of observations")
  }

  # Compute Y-related quantities
  Y_gridlength <- ncol(Y)
  Y_grid_orig <- suppressWarnings(as.numeric(colnames(Y)))
  if (any(is.na(Y_grid_orig)) | is.null(colnames(Y))) {
    Y_grid_orig <- 1:Y_gridlength
  }
  Y_range <- c(0, ((Y_gridlength - 1) / Y_gridlength))
  Y_grid <- seq(Y_range[1], Y_range[2], length.out = Y_gridlength)
  Y_bar <- colMeans(Y)

  # Function to compute ED estimator
  ED_estimator <- function(X_mat) {
    # Compute grid-related quantities
    X_gridlength <- ncol(X_mat)
    X_grid_orig <- suppressWarnings(as.numeric(colnames(X_mat)))
    if (any(is.na(X_grid_orig)) | is.null(colnames(X_mat))) {
      X_grid_orig <- 1:X_gridlength
    }
    X_range <- c(0, ((X_gridlength - 1) / X_gridlength))
    X_grid <- seq(X_range[1], X_range[2], length.out = X_gridlength)
    X_bar <- colMeans(X_mat)

    # Compute covariance and eigendecomposition
    C_hat <- t(X_mat - pracma::repmat(X_bar, T, 1)) %*%
      (Y - pracma::repmat(Y_bar, T, 1)) / T
    D_hat <- C_hat %*% t(C_hat) / Y_gridlength

    # Eigendecomposition
    eigendecomp <- eigen(D_hat)
    psi_hat <- as.matrix(eigendecomp$vectors)
    lambda_hat <- eigendecomp$values

    # Normalize eigenvectors
    norm <- sqrt(colMeans(psi_hat^2))
    psi_hat <- psi_hat / matrix(
      rep(norm, each = nrow(psi_hat)),
      nrow(psi_hat), ncol(psi_hat)
    )
    lambda_hat <- lambda_hat * norm^2

    gscale_X <- sqrt(mean((X_mat - pracma::repmat(X_bar, T, 1))^2))
    gscale_Y <- sqrt(mean((Y - pracma::repmat(Y_bar, T, 1))^2))
    g <- c(1, (2 / pi * atan(gamma * log(T) * lambda_hat / (gscale_X * gscale_Y))))
    K_hat <- which.max(g[-length(g)] - g[-1]) - 1

    list(g = g, K = K_hat)
  }

  # Process functional predictors
  ED_results <- mapply(ED_estimator,
    X = X_functional,
    SIMPLIFY = FALSE
  )

  # Extract K values and name them according to regressors
  K_values <- sapply(functional_vars, function(var) ED_results[[var]]$K)
  names(K_values) <- functional_vars

  # Create plot if requested
  if (plot) {
    for (var in functional_vars) {
      g_values <- ED_results[[var]]$g
      plot(0:(length(g_values) - 1), g_values,
        type = "b",
        xlab = "k", ylab = "g(k)",
        main = paste("g values for", var),
        ylim = c(min(g_values), 1)
      )
      abline(v = K_values[var], col = "red", lty = 2)
      # Pause for user input before showing next plot (except for last plot)
      if (var != tail(functional_vars, 1)) {
        readline(prompt = "Press [Enter] to see next plot")
      }
    }
  }

  # Return only gamma and K values
  return(list(
    gamma = gamma,
    K = K_values
  ))
}
#' Tune gamma hyperparameter in functional eigenvalue difference estimator
#'
#' @description
#' Uses expanding window time-series cross-validation to find the optimal gamma hyperparameter for the factor estimation
#' in functional linear models with \code{\link{fed}}.
#'
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model. Uses the same
#'        format as in \code{\link{flm}}.
#' @param data A list containing the variables in the model, following the same format as in \code{\link{flm}}.
#' @param train.share A number between 0 and 1 specifying the initial proportion of data to use for training.
#'        Default is 0.6 (60% of data for initial training).
#' @param gamma_grid A numeric vector of candidate gamma values to evaluate. Default is 1 to 100 by 1.
#' @param verbose Logical. If TRUE (default), show progress information during computation.
#'
#'
#' @details
#' This function finds the optimal gamma parameter for the \code{\link{fed}} function through
#' cross-validation with an expanding window approach:
#'
#' \enumerate{
#'   \item For each gamma value:
#'     \itemize{
#'       \item Estimate the number of factors (K) for each functional predictor using \code{\link{fed}}
#'       \item Create a formula with only the predictors that have K > 0
#'       \item Initialize training with the first \code{train.share} proportion of observations
#'       \item For each subsequent time point:
#'         \itemize{
#'           \item Fit a model using all data up to the current time point
#'           \item Predict the next observation
#'           \item Compute prediction error
#'         }
#'       \item Average the prediction errors to get the mean squared error (MSE)
#'     }
#'   \item Select the gamma value that minimizes MSE
#' }
#'
#' @return A list containing:
#'   \item{gamma}{The optimal gamma value}
#'   \item{K}{A named vector of the optimal number of factors for each remaining functional predictor}
#'   \item{optimal_formula}{A formula object from the optimal gamma run containing only the predictors with K > 0}
#'   \item{mse}{A vector of mean squared errors for each gamma value}
#'   \item{gamma_grid}{The vector of gamma values evaluated}
#'   \item{removed_vars}{Character vector of variables removed due to having zero factors during the optimal gamma run}
#'
#' @seealso \code{\link{fed}} for factor estimation, \code{\link{flm}} for fitting functional linear models
#'
#' @examples
#' \dontrun{
#' # Simulate functional data
#'
#' # Define a Fourier basis function
#' fourier.basis <- function(evalgrid, nbasis) {
#'   basis <- matrix(NA, length(evalgrid), nbasis)
#'   for (k in 1:nbasis) {
#'     if (k == 1) {
#'       basis[, k] <- rep(1, length(evalgrid))
#'     } else {
#'       if (k %% 2 == 0) {
#'         basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * evalgrid)
#'       } else {
#'         basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * evalgrid)
#'       }
#'     }
#'   }
#'   return(basis)
#' }
#'
#' # Function to generate random bivariate beta coefficients
#' rand.biv.beta <- function(Y.grid, X.grid, K, smoothness_penalty = 1, seed) {
#'   Y.basis <- matrix(NA, length(Y.grid), K)
#'   X.basis <- matrix(NA, length(X.grid), K)
#'   for (k in 1:K) {
#'     if (k == 1) {
#'       Y.basis[, k] <- rep(1, length(Y.grid))
#'       X.basis[, k] <- rep(1, length(X.grid))
#'     } else {
#'       if (k %% 2 == 0) {
#'         Y.basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * Y.grid)
#'         X.basis[, k] <- sqrt(2) * sin(2 * (k / 2) * pi * X.grid)
#'       } else {
#'         Y.basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * Y.grid)
#'         X.basis[, k] <- sqrt(2) * cos(2 * floor(k / 2) * pi * X.grid)
#'       }
#'     }
#'   }
#'
#'   set.seed(seed)
#'
#'   weights <- exp(-seq(0, 2, length.out = K))
#'
#'   # Create coefficient matrix with rank K
#'   A <- matrix(rnorm(K * K), K, K)
#'
#'   # Ensure coefficient matrix has exactly rank K through SVD
#'   svd_A <- svd(A)
#'   A_rank_K <- svd_A$u %*% diag(weights) %*% t(svd_A$v)
#'
#'   # Apply smoothness penalty
#'   A_smooth <- A_rank_K * smoothness_penalty
#'
#'   biv.func <- Y.basis %*% A_smooth %*% t(X.basis)
#'
#'   return(biv.func)
#' }
#'
#' # Set up simulation parameters
#' Y.gridlength <- 200
#' X.gridlength <- 200
#' Y.range <- c(0, ((Y.gridlength - 1) / Y.gridlength))
#' X.range <- c(0, ((X.gridlength - 1) / X.gridlength))
#' Y.grid <- seq(Y.range[1], Y.range[2], length.out = Y.gridlength)
#' X.grid <- seq(X.range[1], X.range[2], length.out = X.gridlength)
#'
#' K <- 3
#' T <- 100 # Number of observations
#'
#' # Generate data
#' four.basis <- fourier.basis(X.grid, 3 * K)
#'
#' f <- matrix(rnorm(n = T * K, 0, 1), T, K)
#' z <- cbind(rep(1, nrow(f)), f)
#' psi <- four.basis[, 1:K]
#' pred.func <- f %*% t(psi)
#' eps <- matrix(rnorm(n = T * 2 * K, 0, 1), T, 2 * K)
#' eps.func <- eps %*% t(four.basis[, (K + 1):(3 * K)])
#'
#' # Create functional predictor
#' X <- pred.func + eps.func
#'
#' # Generate coefficient surface
#' beta <- rand.biv.beta(Y.grid, X.grid, K, seed = 1)
#'
#' intercept <- rep(3, Y.gridlength)
#'
#' # Generate error term
#' u <- matrix(rnorm(n = T * 2 * K, 0, 1), T, 2 * K)
#' u.func <- u %*% t(fourier.basis(Y.grid, 2 * K))
#'
#' # Generate response
#' Y <- intercept + (X %*% beta / X.gridlength) + u.func
#'
#' # Find optimal gamma
#' data <- list(Y = Y, X = X)
#' result <- tune.fed(Y ~ X, data, gamma_grid = seq(1, 10, by = 1))
#'
#' # View results
#' print(result$gamma) # Optimal gamma
#' print(result$K) # Optimal K values
#'
#' # Use optimal parameters in flm
#' model <- flm(result$optimal_formula, data, K = result$K)
#' }
#'
#' @references
#' Otto, S., & Winter, L. (2025). Functional factor regression with an application to electricity price curve modeling.
#' Wu, J. (2018). Eigenvalue difference test for the number of common factors in the approximate factor models. Economics Letters, 169, 63–67.
#'
#' @rdname tune.fed
#' @importFrom stats as.formula
#' @importFrom progress progress_bar
#' @export
tune.fed <- function(formula, data, train.share = 0.6, gamma_grid = seq(1, 100, by = 1), verbose = TRUE) {
  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.list(data)) {
    stop("'data' must be a list containing observations")
  }
  if (!is.numeric(train.share) || train.share <= 0 || train.share >= 1) {
    stop("'train.share' must be a number larger than 0 and smaller than 1")
  }
  if (!is.numeric(gamma_grid) || length(gamma_grid) < 1) {
    stop("'gamma_grid' must be a numeric vector")
  }

  # Extract response variable
  response_var <- all.vars(formula[[2]])
  Y <- if (is.list(data[[response_var]]) && !is.null(data[[response_var]]$type)) {
    data[[response_var]]$data
  } else {
    data[[response_var]]
  }

  # Get total number of observations
  T <- if (is.matrix(Y)) nrow(Y) else length(Y)

  # Calculate initial training size
  train_size <- floor(T * train.share)

  # Initialize storage for results
  mse_results <- numeric(length(gamma_grid))
  K_vectors <- vector("list", length(gamma_grid))
  formulas <- vector("list", length(gamma_grid))
  previous_K <- NULL
  previous_mse <- NULL

  # Initial message
  if (verbose) {
    cat("Starting gamma parameter tuning with", length(gamma_grid), "values...\n")
  }

  # Check for progress package
  has_progress_pkg <- requireNamespace("progress", quietly = TRUE)

  # Create progress bar if available and verbose is TRUE
  if (verbose && has_progress_pkg) {
    pb <- progress::progress_bar$new(
      format = "  Processing [:bar] :percent estimated time: :eta ",
      total = length(gamma_grid),
      clear = FALSE,
      width = 60
    )
  }

  # Parse formula
  form_parts <- as.character(formula)
  response_var <- form_parts[2]
  rhs <- form_parts[3]
  if (rhs == ".") {
    orig_predictors <- predictors <- setdiff(names(data), response_var)
  } else {
    orig_predictors <- all.vars(formula)[-1] # exclude response
  }

  # Helper function to identify scalar variables
  is_scalar <- function(var_data) {
    if (is.list(var_data)) {
      # Check if it's a data list object with explicit type
      if (!is.null(var_data$type)) {
        return(var_data$type == "scalar")
      }
    }
    # Check if it's a data frame with a single column
    if (ncol(as.data.frame(var_data)) == 1) {
      return(TRUE)
    }
    # Default to functional if no type specified and not a single-column list
    return(FALSE)
  }

  # Function to create updated formula
  create_updated_formula <- function(K, orig_predictors, data_list) {
    # Identify scalar and functional predictors
    scalar_vars <- orig_predictors[sapply(data_list[orig_predictors], is_scalar)]
    functional_vars <- setdiff(orig_predictors, scalar_vars)

    # Convert K list to numeric vector if needed
    K_values <- unlist(K)

    # Match K values with functional predictors
    if (length(K_values) != length(functional_vars)) {
      stop("Number of K values does not match number of functional predictors")
    }
    names(K_values) <- functional_vars

    # Get functional terms with non-zero K
    valid_functional_terms <- names(K_values)[K_values > 0]

    if (length(valid_functional_terms) == 0 && length(scalar_vars) == 0) {
      stop("All K values are zero and no scalar predictors remain")
    }

    # Combine scalar variables and valid functional variables for the formula
    all_terms <- c(scalar_vars, valid_functional_terms)

    # Create new formula
    as.formula(paste(
      response_var,
      "~",
      paste(all_terms, collapse = " + ")
    ))
  }

  # Loop through gamma grid
  for (g_idx in seq_along(gamma_grid)) {
    gamma <- gamma_grid[g_idx]

    # Update progress based on method available
    if (verbose) {
      if (has_progress_pkg) {
        pb$tick()
      } else {
        cat(sprintf("Processing gamma = %.3f (%d/%d)\n", gamma, g_idx, length(gamma_grid)))
      }
    }

    # Run FED test
    fed_result <- fed(formula, data, gamma = gamma, plot = FALSE)

    # Store K vector
    current_K <- fed_result$K[fed_result$K > 0]
    K_vectors[[g_idx]] <- current_K

    # Create updated formula based on non-zero K values
    current_formula <- try(create_updated_formula(fed_result$K, orig_predictors, data), silent = TRUE)

    if (inherits(current_formula, "try-error")) {
      if (verbose) {
        cat("Error creating formula for gamma =", gamma, ". Skipping.\n")
      }
      mse_results[g_idx] <- Inf
      next
    }

    formulas[[g_idx]] <- current_formula

    # Get predictors in current formula
    current_predictors <- all.vars(current_formula)[-1]

    # Check if K vector is different from previous iteration
    if (g_idx > 1 && identical(current_K, previous_K)) {
      # If K is the same, use previous MSE
      if (verbose && !has_progress_pkg) {
        cat("  Same K vector as previous gamma, reusing MSE\n")
      }
      mse_results[g_idx] <- previous_mse
    } else {
      # Initialize prediction errors matrix
      pred_errors <- matrix(0, nrow = T - train_size, ncol = ncol(Y))

      # Expanding window validation
      for (t in train_size:(T - 1)) {
        # Create training subset
        train_data <- lapply(data, function(x) {
          if (is.matrix(x)) {
            x[1:t, , drop = FALSE]
          } else if (is.list(x) && !is.null(x$type) && x$type == "scalar") {
            list(type = "scalar", data = x$data[1:t])
          } else if (is.numeric(x)) {
            x[1:t]
          } else {
            x
          }
        })

        # Fit model
        model <- flm(current_formula, train_data, K = current_K, inference = FALSE)

        # Prepare validation data
        valid_data <- lapply(current_predictors, function(var) {
          x <- data[[var]]
          if (is.matrix(x)) {
            x[(t + 1):(t + 1), , drop = FALSE]
          } else if (is.list(x) && !is.null(x$type) && x$type == "scalar") {
            list(type = "scalar", data = x$data[t + 1])
          } else if (is.numeric(x)) {
            x[t + 1]
          } else {
            x
          }
        })
        names(valid_data) <- current_predictors

        # Make prediction
        pred <- predict(model, valid_data)
        actual <- if (is.matrix(Y)) Y[t + 1, ] else Y[t + 1]
        pred_errors[t - train_size + 1, ] <- (actual - pred)^2
      }

      # Store MSE
      mse_results[g_idx] <- mean(pred_errors)

      # Store current values for next iteration
      previous_mse <- mse_results[g_idx]
    }

    previous_K <- current_K
  }

  # Find best gamma (minimum MSE)
  best_idx <- which.min(mse_results)

  # Get final results
  best_formula <- formulas[[best_idx]]
  best_K <- K_vectors[[best_idx]]
  removed_vars <- setdiff(orig_predictors, all.vars(best_formula)[-1])

  if (verbose) {
    cat("\nTuning complete:\n")
    cat(sprintf("  Optimal gamma: %.3f\n", gamma_grid[best_idx]))
    cat(sprintf("  Optimal K values: %s\n", paste(best_K, collapse = ", ")))
    cat(sprintf("  MSE: %.6f\n", mse_results[best_idx]))

    if (length(removed_vars) > 0) {
      cat("  Variables removed due to zero factors:", paste(removed_vars, collapse = ", "), "\n")
    }
  }

  # Return results
  list(
    gamma = gamma_grid[best_idx],
    K = best_K,
    optimal_formula = best_formula,
    mse = mse_results,
    gamma_grid = gamma_grid,
    removed_vars = removed_vars
  )
}
