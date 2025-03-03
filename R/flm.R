#' Fitting fully functional linear models
#'
#' @description
#' Fits a linear model with functional response and functional/scalar predictors.
#' This implementation is based on the functional factor regression approach which enables statistical
#' inference in function-on-function linear regressions.
#'
#' @param formula an object of class \code{\link[stats]{formula}} which is a symbolic description of the model to be fitted.
#'        A model takes the form \code{response ~ regressors}, where \code{response} is the functional response matrix and \code{regressors}
#'        are a series of functional and/or scalar response matrices/vectors. A "regressors" specification of the form
#'        \code{regressor1 + regressor2} indicates a regression on both variables. A formula specification \code{response ~ .} implies
#'        a regression on all other variables in \code{data} besides \code{response}.
#' @param data an object of class \code{list} which contains the variables in the model. Each named object within the list is considered
#'        either a functional or a scalar variable. Functional variables are matrices with the number of sample observations T equal the
#'        number of rows and the number of functional observations on an equidistant grid equal the number of columns. Scalar variables
#'        are vectors of length equal to the number of sample observations T, or matrices of column length 1. Additionally, one
#'        can define the \code{type} of each object within the list as either \code{functional} or \code{scalar} for improved reproducibility. While
#'        the number of sample observations must be the same for all objects, the number of functional observations might differ from
#'        functional variable to functional variable.
#' @param K a vector of positive integers specifying the number of factors for each functional regressor.
#'        Must be the same length as the number of functional predictors in the model.
#'        Names should match the functional predictor names in \code{data}. For anonymous vectors,
#'        the values are matched according to the order of regressors in the \code{formula}.
#'        See \code{\link{fed}} and \code{\link{tune.fed}} for consistent estimation of the number of factors.
#' @param conf.level confidence level (usually between 0.9 and 0.999) used for calculating the confidence regions for
#'        all bivariate function coefficients. Default is 0.95 (95% confidence).
#' @param inference logical. If \code{FALSE}, only the point estimates of all regression objects are computed which speeds up the runtime
#'        significantly. This might be useful in combination with \code{\link{predict.flm}} for large forecasting loops. Default is \code{TRUE}. Note that
#'        \code{inference = FALSE} leads to limited functionality of \code{\link{plot.summary.flm}}.
#'
#' @details
#' The function implements the functional factor regression approach for statistical inference in function-on-function linear models,
#' introduced by Otto & Winter (2025). It handles both functional and scalar predictors and provides confidence regions when \code{inference = TRUE}.
#' Consult the original paper for details about the method.
#'
#' Formula specification examples:
#' \itemize{
#'   \item \code{Y ~ X} - Simple functional regression with one functional predictor
#'   \item \code{Y ~ X1 + X2} - Multiple functional predictors
#'   \item \code{Y ~ X + w} - Mixed functional and scalar predictors
#'   \item \code{Y ~ .} - Use all variables in \code{data} except the response
#' }
#'
#' @return An object of class "flm" containing:
#'   \item{info}{List with the matched call, number of observations, and grid information about the data:
#'     \itemize{
#'       \item{call}{The matched function call}
#'       \item{n_obs}{Number of observations}
#'       \item{grid_info}{Grid information for the response and predictors}
#'     }
#'   }
#'   \item{coefficients}{List with three components:
#'     \itemize{
#'       \item{beta}{Functional coefficient surfaces for each predictor}
#'       \item{B}{Matrix of coefficients from finite regression representation}
#'       \item{scalar}{Univariate functional coefficients for the intercept and scalar predictors}
#'     }
#'   }
#'   \item{residuals}{List containing:
#'     \itemize{
#'       \item{Y}{The regression residuals}
#'       \item{X}{The regressors' factor model equation errors}
#'     }
#'   }
#'   \item{fitted}{List containing:
#'     \itemize{
#'       \item{Y}{The fitted values of the response}
#'       \item{factors}{Information about the fitted factor models for the regressors}
#'     }
#'   }
#'   \item{inference}{When \code{inference=TRUE}, a list containing:
#'     \itemize{
#'       \item{covariance}{Bivariate covariance functions}
#'       \item{t_values}{t-values for coefficient surfaces}
#'       \item{conf_bands}{Confidence regions for coefficient surfaces}
#'       \item{conf_level}{The specified confidence level}
#'     }
#'   }
#'
#' @section Methods:
#' The following methods are available for "flm" objects:
#' \itemize{
#'   \item \code{\link{print.flm}} for basic model output
#'   \item \code{\link{summary.flm}} for detailed model summaries
#'   \item \code{\link{plot.summary.flm}} for visualizing model components
#'   \item \code{\link{predict.flm}} for generating predictions
#' }
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
#' # Create scalar predictor
#' w <- as.matrix(rnorm(n = T, 0, 1))
#' alpha <- cos(Y.grid)
#'
#' intercept <- rep(3, Y.gridlength)
#'
#' # Generate error term
#' u <- matrix(rnorm(n = T * 2 * K, 0, 1), T, 2 * K)
#' u.func <- u %*% t(fourier.basis(Y.grid, 2 * K))
#'
#' # Generate response
#' Y <- intercept + (w %*% t(alpha)) + (X %*% beta / X.gridlength) + u.func
#'
#' # Fit functional factor regression model
#' data <- list(Y = Y, X = X, w = list(type = "scalar", data = w))
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3), conf.level = 0.99, inference = TRUE)
#'
#' # Examine model results
#' summary(ffr_model)
#' plot(summary(ffr_model))
#'
#' # Make predictions
#' pred <- predict(ffr_model, newdata = data)
#' }
#'
#' @references
#' Otto, S., & Winter, L. (2025). Functional factor regression with an application to electricity price curve modeling.
#'
#' @seealso
#' \code{\link{summary.flm}} for model summaries
#' \code{\link{plot.summary.flm}} for summary plots
#' \code{\link{predict.flm}} for prediction
#' \code{\link{fed}} and \code{\link{tune.fed}} for factor estimation
#' \code{\link[stats]{formula}} for formula specification details
#'
#' @rdname flm
#' @importFrom pracma repmat
#' @keywords models regression functional inference
#' @export
flm <- function(formula, data, K, conf.level = 0.95, inference = TRUE) {
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

  # Validate K for functional predictors only
  if (missing(K)) {
    stop("'K' must be provided to specify the number of factors for functional regressors")
  }
  if (!is.numeric(K) || any(K <= 0) || any(K %% 1 != 0)) {
    stop("'K' must be a vector of positive integers")
  }
  if (length(K) != length(functional_vars)) {
    stop(sprintf(
      "Length of K (%d) must match number of functional predictors (%d)",
      length(K), length(functional_vars)
    ))
  }

  # Name the K vector if it's not already named
  if (is.null(names(K))) {
    names(K) <- functional_vars
  } else {
    if (!all(sort(names(K)) == sort(functional_vars))) {
      stop("Names of K must match functional predictor variables")
    }
  }

  # Extract response and predictor data
  Y <- if (is_scalar(data[[response_var]])) data[[response_var]]$data else data[[response_var]]

  # Extract and validate functional predictors
  X_functional <- lapply(functional_vars, function(var) {
    if (is_scalar(data[[var]])) data[[var]]$data else data[[var]]
  })
  names(X_functional) <- functional_vars

  # Extract scalar predictors
  X_scalar <- lapply(scalar_vars, function(var) {
    var_data <- data[[var]]
    if (is_scalar(var_data)) {
      # Check if it's a data list object with $data structure
      if (is.list(var_data) && !is.null(var_data$data)) {
        return(var_data$data)
      }
      # Otherwise return the data directly
      return(var_data)
    }
    return(var_data)
  })
  names(X_scalar) <- scalar_vars

  # Validate dimensions
  T <- if (is.matrix(Y)) nrow(Y) else length(Y)

  # Validate functional predictor dimensions
  if (!all(sapply(X_functional, nrow) == T)) {
    stop("All functional variables must have the same number of observations")
  }

  # Validate scalar predictor dimensions
  if (!all(sapply(X_scalar, length) == T)) {
    stop("All scalar variables must have the same number of observations")
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

  # Initialize z_hat with intercept
  z_hat <- matrix(1, nrow = T, ncol = 1)

  # Add scalar predictors directly to z_hat
  for (var in scalar_vars) {
    z_hat <- cbind(z_hat, as.matrix(X_scalar[[var]]))
  }

  # Function to process each functional regressor
  process_regressor <- function(X_mat, K) {
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
    psi_hat <- as.matrix(eigendecomp$vectors[, 1:K])
    lambda_hat <- eigendecomp$values[1:K]

    # Normalize eigenvectors
    norm <- sqrt(colMeans(psi_hat^2))
    psi_hat <- psi_hat / matrix(
      rep(norm, each = nrow(psi_hat)),
      nrow(psi_hat), ncol(psi_hat)
    )
    lambda_hat <- lambda_hat * norm^2

    # Compute scores
    f_hat <- (X_mat - pracma::repmat(X_bar, T, 1)) %*% psi_hat / X_gridlength

    # Compute regressor errors
    epsilon_hat <- X_mat - pracma::repmat(X_bar, T, 1) - f_hat %*% t(psi_hat)

    # Return all computed quantities
    list(
      data = X_mat,
      grid = X_grid_orig,
      mean = X_bar,
      psi = psi_hat,
      lambda = lambda_hat,
      f = f_hat,
      epsilon = epsilon_hat
    )
  }

  # Process functional predictors
  processed_X <- mapply(process_regressor,
    X = X_functional,
    K = K[functional_vars],
    SIMPLIFY = FALSE
  )

  # Consolidate processed matrices
  Psi_hat <- lapply(processed_X, function(x) x$psi)
  Lambda_hat <- lapply(processed_X, function(x) x$lambda)
  F_hat <- lapply(processed_X, function(x) x$f)
  epsilon_hat <- lapply(processed_X, function(x) x$epsilon)

  # Add functional components to z_hat
  for (f in F_hat) {
    z_hat <- cbind(z_hat, f)
  }

  # Modify get_subset to account for scalar predictors
  get_subset <- function(j, K_values) {
    J <- length(K_values)
    N <- 2 + length(scalar_vars)
    if (j < 1 || j > J) {
      stop(sprintf("j must be between 1 and %d", J))
    }
    start_idx <- N + if (j > 1) sum(K_values[1:(j - 1)]) else 0
    K_j <- K_values[j]
    rows <- start_idx:(start_idx + K_j - 1)
    return(rows)
  }

  # Estimate regression coefficients
  B_hat <- solve(t(z_hat) %*% z_hat) %*% t(z_hat) %*% Y

  # Calculate fitted values and residuals
  Y_hat <- z_hat %*% B_hat
  u_hat <- Y - Y_hat

  # Recover bivariate functional regression coefficients
  beta_hat <- list()
  for (j in seq_along(functional_vars)) {
    subset_indices <- get_subset(j, K)
    beta_hat[[functional_vars[j]]] <- t(B_hat[subset_indices, , drop = F]) %*% t(Psi_hat[[j]])
  }

  # Only compute inference-related quantities if inference = TRUE
  inference_results <- if (inference) {
    # Compute covariance-related quantities
    Q_hat <- t(z_hat) %*% z_hat / T
    Q_hat_inv <- solve(Q_hat)
    z_bar <- as.matrix(colMeans(z_hat))

    # Process each regressor's covariance components
    process_covariance_components <- function(j, F_hat_j, Psi_hat_j, Lambda_hat_j, epsilon_hat_j) {
      zF_bar_j <- t(z_hat) %*% F_hat_j / T
      gamma_hat_j <- t(F_hat_j) %*% (Y - pracma::repmat(Y_bar, T, 1)) / T
      Fy_bar_j <- t(F_hat_j) %*% t(gamma_hat_j %*% t(Y - pracma::repmat(Y_bar, T, 1)) / Y_gridlength) / T

      list(
        zF_bar = zF_bar_j,
        gamma_hat = gamma_hat_j,
        Fy_bar = Fy_bar_j
      )
    }

    # Apply covariance processing to all regressors
    cov_components <- mapply(process_covariance_components,
      j = seq_along(functional_vars),
      F_hat_j = F_hat,
      Psi_hat_j = Psi_hat,
      Lambda_hat_j = Lambda_hat,
      epsilon_hat_j = epsilon_hat,
      SIMPLIFY = FALSE
    )

    zF_bar <- lapply(cov_components, function(x) x$zF_bar)
    gamma_hat <- lapply(cov_components, function(x) x$gamma_hat)
    Fy_bar <- lapply(cov_components, function(x) x$Fy_bar)

    # Function to compute G matrices for a specific regressor and time point
    compute_G_matrix <- function(k, F_hat_k, gamma_hat_k, Lambda_hat_k, Fy_bar_k, t) {
      K_k <- length(Lambda_hat_k)
      G_hat <- matrix(0, K_k, K_k)

      for (l in 1:K_k) {
        for (m in 1:K_k) {
          if (l != m) {
            G_hat[l, m] <- (F_hat_k[t, m] * mean(gamma_hat_k[l, ] * (Y[t, ] - Y_bar)) - Fy_bar_k[m, l] +
              F_hat_k[t, l] * mean(gamma_hat_k[m, ] * (Y[t, ] - Y_bar)) - Fy_bar_k[l, m]) /
              (Lambda_hat_k[m] - Lambda_hat_k[l])
          }
        }
      }
      return(G_hat)
    }

    # Compute Omega for each regressor
    compute_omega <- function(j) {
      Omega_j <- matrix(0, Y_gridlength, ncol(X_functional[[functional_vars[j]]]))
      Qinv_subset <- get_subset(j, K)
      Q_hat_inv_sub <- Q_hat_inv[Qinv_subset, ]

      for (t in 1:T) {
        # Compute A1
        A1_t <- t(Q_hat_inv_sub %*% t(z_hat[t, , drop = FALSE]) %*% u_hat[t, , drop = FALSE]) %*%
          t(Psi_hat[[j]])

        # Compute A2
        temp_A2_t <- matrix(0, nrow(z_bar), Y_gridlength)
        for (k in seq_along(functional_vars)) {
          G_hat_k <- compute_G_matrix(
            k, F_hat[[k]], gamma_hat[[k]],
            Lambda_hat[[k]], Fy_bar[[k]], t
          )
          Bhatk_subset <- get_subset(k, K)
          temp_A2_t <- temp_A2_t +
            (z_bar %*% F_hat[[k]][t, , drop = FALSE] - zF_bar[[k]] %*% G_hat_k) %*%
            B_hat[Bhatk_subset, ]

          if (k == j) {
            G_hat_j <- G_hat_k
          }
        }
        A2_t <- Psi_hat[[j]] %*% Q_hat_inv_sub %*% temp_A2_t

        # Compute A3
        Bhatj_subset <- get_subset(j, K)
        A3_t <- Psi_hat[[j]] %*% G_hat_j %*% B_hat[Bhatj_subset, ]

        # Compute A4
        A4_t <- t(epsilon_hat[[j]][t, , drop = FALSE]) %*%
          t((gamma_hat[[j]] %*% t(Y[t, , drop = FALSE] - Y_bar) / Y_gridlength) /
            Lambda_hat[[j]]) %*% B_hat[Bhatj_subset, ]

        # Combine all components
        omega_t <- t(A1_t) + A2_t + A3_t + A4_t
        Omega_j <- Omega_j + (t(omega_t^2) / T)
      }
      return(Omega_j)
    }

    # Compute Omega for all regressors
    Omega <- lapply(seq_along(functional_vars), compute_omega)
    names(Omega) <- functional_vars

    # Compute t-values for each regressor
    t_values <- mapply(function(beta_j, omega_j) {
      sqrt(T) * beta_j / sqrt(omega_j)
    }, beta_j = beta_hat, omega_j = Omega, SIMPLIFY = FALSE)
    names(t_values) <- functional_vars

    # Compute confidence bands for each regressor
    confidence_bands <- mapply(function(beta_j, omega_j) {
      z_value <- qnorm(1 - ((1 - conf.level) / 2))
      std_error <- sqrt(omega_j / T)
      list(
        lower = beta_j - z_value * std_error,
        upper = beta_j + z_value * std_error
      )
    }, beta_j = beta_hat, omega_j = Omega, SIMPLIFY = FALSE)
    names(confidence_bands) <- functional_vars

    # Return inference results
    list(
      covariance = Omega,
      t_values = t_values,
      conf_bands = confidence_bands,
      conf_level = conf.level
    )
  } else {
    # Return NULL if inference is FALSE
    NULL
  }

  # Organize input information into a single list
  model_info <- list(
    call = match.call(),
    n_obs = T,
    grid_info = list(
      Y = list(
        data = Y,
        grid = Y_grid_orig,
        length = Y_gridlength,
        mean = Y_bar
      ),
      X = c(
        # Scalar predictors info
        setNames(
          lapply(scalar_vars, function(var) {
            list(
              name = var,
              type = "scalar"
            )
          }),
          scalar_vars
        ),
        # Functional predictors info
        setNames(
          lapply(seq_along(functional_vars), function(j) {
            list(
              name = functional_vars[j],
              type = "functional",
              data = processed_X[[j]]$data,
              grid = processed_X[[j]]$grid,
              length = ncol(X_functional[[j]]),
              mean = processed_X[[j]]$mean,
              K = K[j]
            )
          }),
          functional_vars
        )
      )
    )
  )

  # Create final model object
  model <- list(
    info = model_info,
    coefficients = list(
      beta = beta_hat, # functional coefficients
      B = B_hat, # basis coefficients
      scalar = if (length(scalar_vars) > 0) {
        # Extract scalar coefficients (intercept and scalar predictors)
        scalar_coef <- B_hat[1:(1 + length(scalar_vars)), , drop = FALSE]
        rownames(scalar_coef) <- c("intercept", scalar_vars)
        scalar_coef
      } else {
        # Only intercept if no scalar predictors
        structure(B_hat[1, , drop = FALSE],
          dimnames = list("intercept", NULL)
        )
      }
    ),
    residuals = list(
      Y = u_hat,
      X = lapply(processed_X, function(x) x$epsilon)
    ),
    fitted = list(
      Y = Y_hat,
      factors = list(
        Psi = lapply(processed_X, function(x) x$psi),
        Lambda = lapply(processed_X, function(x) x$lambda),
        F = F_hat
      )
    ),
    inference = inference_results # Add inference results (NULL if inference = FALSE)
  )

  class(model) <- c("flm", "regression")
  return(model)
}
