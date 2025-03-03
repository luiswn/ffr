#' Natural spline preprocessing of functional data
#'
#' @description
#' Transforms discrete point observations of functional data into functional data evaluated on a dense equidistant working grid.
#' The function uses a direct interpolation approach, using natural splines to interpolate missing values
#'
#' @param data A numeric matrix or a list of numeric matrices. For matrices, rows represent observations
#'        and columns represent points on the functional grid. Scalar variables (single-column matrices)
#'        are preserved unchanged.
#' @param grid_length An integer specifying the desired working grid length for the functional data.
#'
#' @details
#' This function transforms discretely observed functional data objects to densely evaluated variables.
#' For each observation (row):
#'
#' \enumerate{
#'   \item A natural spline is fitted to the original data points
#'   \item The fitted spline is evaluated at the new grid points
#'   \item Missing values (NA) in the original data are handled appropriately
#' }
#'
#' The function preserves scalar variables (single-column matrices) without modification and
#' maintains the structure of the input (matrix or list).
#'
#' @return
#' If \code{data} is a matrix, returns a matrix with the same number of rows and \code{grid_length} columns.
#' If \code{data} is a list, returns a list with the same structure, where each functional matrix has been
#' resampled to have \code{grid_length} columns. Column names of the output matrices represent the new grid points.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' X1 <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50) # 100 observations, 50 grid points
#' X2 <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30) # 100 observations, 30 grid points
#' scalar_var <- matrix(rnorm(100), nrow = 100, ncol = 1) # scalar variable
#'
#' # Process a single matrix
#' X1_processed <- fpreprocess(X1, grid_length = 100)
#'
#' # Process a list of matrices
#' data_list <- list(X1 = X1, X2 = X2, scalar = scalar_var)
#' processed_data <- fpreprocess(data_list, grid_length = 100)
#' }
#'
#' @rdname fpreprocess
#' @importFrom splines ns
#' @importFrom stats lm
#' @export
#'
fpreprocess <- function(data, grid_length) {
  # Helper function to check if a variable is scalar (single column)
  is_scalar <- function(mat) {
    return(is.matrix(mat) && ncol(mat) == 1)
  }

  # Helper function to process a single matrix
  process_matrix <- function(mat, new_length) {
    if (!is.numeric(mat) || !is.matrix(mat)) {
      stop("Each element must be a numeric matrix")
    }

    # If scalar (single column), return unchanged
    if (is_scalar(mat)) {
      return(mat)
    }

    # Create new working grid
    old_grid <- seq_len(ncol(mat))
    new_grid <- seq(1, ncol(mat), length.out = new_length)

    # Initialize output matrix
    splinefit <- matrix(nrow = nrow(mat), ncol = new_length)

    # Fit natural splines for each row
    splinefit <- t(apply(mat, 1, function(row) {
      valid_idx <- !is.na(row)
      thisobsgrid <- old_grid[valid_idx]
      thisdata <- row[valid_idx]

      # Fit spline
      splinebasis <- splines::ns(thisobsgrid,
        knots = thisobsgrid[-length(thisobsgrid)],
        Boundary.knots = c(thisobsgrid[1], thisobsgrid[length(thisobsgrid)])
      )

      coef <- lm(thisdata ~ splinebasis - 1)$coefficients

      # Interpolate on new grid
      densebasis <- splines::ns(new_grid,
        knots = thisobsgrid[-length(thisobsgrid)],
        Boundary.knots = c(thisobsgrid[1], thisobsgrid[length(thisobsgrid)])
      )

      densebasis %*% coef
    }))

    # Preserve row names if they exist
    if (!is.null(rownames(mat))) {
      rownames(splinefit) <- rownames(mat)
    }

    # Create column names for new grid
    colnames(splinefit) <- new_grid

    return(splinefit)
  }

  # Main function logic
  if (is.list(data) && !is.matrix(data) && !is.data.frame(data)) {
    # Process list of matrices
    result <- lapply(data, function(mat) {
      process_matrix(mat, grid_length)
    })
    # Preserve list names
    names(result) <- names(data)
  } else {
    # Process single matrix
    result <- process_matrix(data, grid_length)
  }

  return(result)
}
