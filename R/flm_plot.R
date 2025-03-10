#' Plot method for summary.flm objects
#'
#' @description
#' Creates visualizations of functional linear model components including coefficient surfaces,
#' t-values, p-values, and R-squared goodness-of-fit measure.
#'
#' @param x An object of class "summary.flm", typically the result of a call to \code{\link{summary.flm}}.
#' @param predictor Character string specifying which predictor to visualize. Use "intercept" for the
#'        intercept function. Not required when \code{which = "R2"}.
#' @param which Character string specifying the type of plot to create. Options are:
#'        \itemize{
#'          \item{\code{"beta"}}{Coefficient surface (default)}
#'          \item{\code{"t"}}{t-values surface}
#'          \item{\code{"p"}}{p-values surface with significance contours}
#'          \item{\code{"beta_3D"}}{3D visualization of coefficient surface (requires plotly)}
#'          \item{\code{"R2"}}{Functional R-squared across the response grid}
#'        }
#' @param conf.region Logical. If \code{TRUE}, confidence regions are added to coefficient plots when available.
#'        Default is \code{FALSE}.
#' @param ... Additional arguments passed to summary plotting functions.
#'
#' @details
#' This function provides a variety of visualization options for interpreting functional linear models:
#'
#' \subsection{Plot Types}{
#'   \itemize{
#'     \item{\strong{Coefficient plots} (\code{which = "beta"})}{
#'       For functional predictors: Creates a heatmap of the bivariate coefficient surface.
#'       For scalar predictors: Creates a line plot of the coefficient function.
#'       For intercept: Creates a line plot of the intercept function.
#'     }
#'     \item{\strong{3D coefficient plots} (\code{which = "beta_3D"})}{
#'       Creates an interactive 3D surface plot of the coefficient (requires plotly package).
#'       Only available for functional predictors.
#'     }
#'     \item{\strong{t-value plots} (\code{which = "t"})}{
#'       Visualizes the pointwise t-statistic values across the grid of the coefficient surface.
#'       Only available for functional predictors.
#'     }
#'     \item{\strong{p-value plots} (\code{which = "p"})}{
#'       Displays pointwise p-values with contour lines at significance levels 0.01, 0.05, and 0.1.
#'       Only available for functional predictors.
#'     }
#'     \item{\strong{R-squared plot} (\code{which = "R2"})}{
#'       Shows how the R-squared measure varies across the response grid, indicating where the model
#'       fits better or worse.
#'     }
#'   }
#' }
#'
#' \subsection{Confidence Regions}{
#'   When \code{conf.region = TRUE} and the model was fitted with \code{inference = TRUE},
#'   confidence bands will be displayed:
#'   \itemize{
#'     \item For 2D coefficient plots: Three panels showing the point estimate, lower bound, and upper bound.
#'     \item For 3D plots: Additional semi-transparent surfaces for the confidence regions.
#'   }
#' }
#'
#' \subsection{Required Packages}{
#'   Some plot types require additional packages:
#'   \itemize{
#'     \item The \code{fields} package is required for heatmap plots.
#'     \item The \code{scales} package is required for color scales.
#'     \item The \code{plotly} package is required for 3D plots.
#'   }
#' }
#'
#' @return
#' For \code{which = "beta_3D"}, returns a plotly object. For all other plot types,
#' returns the plot invisibly (called for creating a plot).
#'
#' @seealso \code{\link{flm}}, \code{\link{summary.flm}}
#'
#' @examples
#' \dontrun{
#' # Fit a functional linear model
#' ffr_model <- flm(Y ~ X + w, data = data, K = c(3), inference = TRUE)
#' model_summary <- summary(ffr_model)
#'
#' # Plot coefficient surface for a functional predictor
#' plot(model_summary, predictor = "X", which = "beta")
#'
#' # Plot with confidence regions
#' plot(model_summary, predictor = "X", which = "beta", conf.region = TRUE)
#'
#' # Plot 3D visualization of coefficient
#' plot(model_summary, predictor = "X", which = "beta_3D")
#'
#' # Plot statistical significance
#' plot(model_summary, predictor = "X", which = "t") # t-values
#' plot(model_summary, predictor = "X", which = "p") # p-values
#'
#' # Plot coefficient function for a scalar predictor
#' plot(model_summary, predictor = "w", which = "beta")
#'
#' # Plot intercept function
#' plot(model_summary, predictor = "intercept", which = "beta")
#'
#' # Plot R-squared across response grid
#' plot(model_summary, which = "R2")
#' }
#'
#' @rdname plot.summary.flm
#' @method plot summary.flm
#' @importFrom magrittr %>%
#' @importFrom plotly plot_ly add_surface layout
#' @importFrom fields image.plot imageplot.setup
#' @importFrom scales alpha
#' @importFrom graphics abline contour image legend matplot par
#' @importFrom stats pt
#' @importFrom grDevices colorRampPalette
#' @export
plot.summary.flm <- function(x, predictor = NULL,
                             which = c("beta", "t", "p", "beta_3D", "R2"),
                             conf.region = FALSE, ...) {
  # Input validation
  which <- match.arg(which)

  # Handle R2 case separately and early
  if (which == "R2") {
    matplot(x$info$grid_info$Y$grid, x$stats$R_squared_func,
      type = "l", lty = 1,
      xlab = paste("Grid points", x$info$call$formula[[2]]), ylab = "",
      main = "Functional R-squared goodness-of-fit measure"
    )
    return(invisible())
  }

  # For all other plots, we need a valid predictor
  if (is.null(predictor)) {
    stop("Predictor must be specified for plots other than R2")
  }

  # Special handling for intercept
  if (predictor == "intercept" || predictor == "(Intercept)") {
    if (which != "beta") {
      warning("Only 'beta' plotting is available for intercept")
      return(NULL)
    }

    # Create plot for intercept across Y grid
    matplot(x$info$grid_info$Y$grid, as.vector(x$coefficients$intercept),
      type = "l", lty = 1,
      xlab = paste("Grid points", x$info$call$formula[[2]]), ylab = "Coefficient value",
      main = "Intercept estimate"
    )
    abline(h = 0, lty = 2, col = "gray")

    return(invisible())
  }

  if (!(predictor %in% names(x$info$grid_info$X))) {
    stop("Predictor not found in model")
  }

  # Check if inference results are available when conf.region is requested
  if (conf.region && is.null(x$inference)) {
    warning("Confidence regions cannot be plotted as the model was fitted with inference = FALSE")
    conf.region <- FALSE
  }

  # Check if predictor is scalar
  is_scalar <- x$info$grid_info$X[[predictor]]$type == "scalar"

  # Extract grids
  Y_grid <- x$info$grid_info$Y$grid
  X_grid <- x$info$grid_info$X[[predictor]]$grid

  if (is_scalar) {
    if (which != "beta") {
      warning("Only 'beta' plotting is available for scalar predictors")
      return(NULL)
    }

    # Plot scalar coefficient
    coef_idx <- which(rownames(x$coefficients$scalar) == predictor)
    coef_value <- as.vector(x$coefficients$scalar[coef_idx, ])

    # Create plot for scalar coefficient across Y grid
    matplot(Y_grid, coef_value,
      type = "l", lty = 1,
      xlab = paste("Grid points", x$info$call$formula[[2]]), ylab = "Coefficient value",
      main = paste("Coefficient estimate for", predictor)
    )
    abline(h = 0, lty = 2, col = "gray")

    return(invisible())
  }

  # For functional predictors, load required packages
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Please install the 'fields' package to use this function")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Please install the 'scales' package to use this function")
  }
  if (which == "beta_3D" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Please install the 'plotly' package to use the 3D plot")
  }


  # Ensure grids are sorted
  if (is.unsorted(X_grid)) {
    oX <- order(X_grid)
    X_grid <- X_grid[oX]
  } else {
    oX <- seq_along(X_grid)
  }
  if (is.unsorted(Y_grid)) {
    oY <- order(Y_grid)
    Y_grid <- Y_grid[oY]
  } else {
    oY <- seq_along(Y_grid)
  }

  # Extract coefficient and confidence bands
  beta <- x$coefficients$beta[[predictor]][oY, oX, drop = FALSE]
  if (conf.region) {
    lower <- x$inference$conf_bands[[predictor]]$lower[oY, oX, drop = FALSE]
    upper <- x$inference$conf_bands[[predictor]]$upper[oY, oX, drop = FALSE]
    zlim_beta <- max(abs(c(beta, lower, upper)), na.rm = TRUE) * c(-1, 1)
  } else {
    zlim_beta <- max(abs(beta), na.rm = TRUE) * c(-1, 1)
  }
  cols_beta <- colorRampPalette(c("blue", "white", "red"))(100)

  if (which == "beta_3D") {
    # Create meshgrid for 3D plotting
    grid_x <- matrix(rep(X_grid, length(Y_grid)), ncol = length(Y_grid))
    grid_y <- t(matrix(rep(Y_grid, length(X_grid)), ncol = length(X_grid)))

    p <- plotly::plot_ly() %>%
      plotly::add_surface(
        x = ~grid_x,
        y = ~grid_y,
        z = ~ t(beta),
        colorscale = list(
          c(0, 0.5, 1),
          c("blue", "white", "red")
        ),
        cmin = zlim_beta[1],
        cmax = zlim_beta[2],
        name = "beta(r,s)"
      )

    if (conf.region) {
      p <- p %>%
        plotly::add_surface(
          x = ~grid_x,
          y = ~grid_y,
          z = ~ t(lower),
          colorscale = list(c(0, 1), c("grey", "grey")),
          opacity = 0.3,
          showscale = FALSE,
          name = paste0("Lower ", x$inference$conf_level * 100, "% Band")
        ) %>%
        plotly::add_surface(
          x = ~grid_x,
          y = ~grid_y,
          z = ~ t(upper),
          colorscale = list(c(0, 1), c("grey", "grey")),
          opacity = 0.3,
          showscale = FALSE,
          name = paste0("Upper ", x$inference$conf_level * 100, "% Band")
        )
    }

    p <- p %>% plotly::layout(
      title = list(
        text = paste0(
          "3D visualization of beta(r,s)",
          if (conf.region) paste0(" with ", x$inference$conf_level * 100, "% confidence bands") else ""
        ),
        x = 0.5,
        y = 0.95
      ),
      scene = list(
        xaxis = list(title = paste("Grid points", predictor)),
        yaxis = list(title = paste("Grid points", x$info$call$formula[[2]])),
        zaxis = list(title = "Value"),
        camera = list(
          eye = list(x = -2, y = -2, z = 2),
          center = list(x = 0.5, y = 0.5, z = 0),
          up = list(x = 0, y = 0, z = 1)
        )
      )
    )

    return(p)
  } else if (which == "beta") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (conf.region) {
      # Set up 3-panel plot
      par(mfrow = c(1, 3), mar = c(5, 4, 4, 4))

      # Create a common color legend
      legend_info <- fields::imageplot.setup(
        zlim = zlim_beta,
        col = cols_beta
      )

      # Plot estimated coefficient
      image(
        x = X_grid, y = Y_grid, z = t(beta),
        zlim = zlim_beta, col = cols_beta,
        xlab = paste("Grid points", predictor),
        ylab = paste("Grid points", x$info$call$formula[[2]]),
        main = expression(paste("Estimated ", beta(r, s)))
      )

      # Plot lower confidence band
      image(
        x = X_grid, y = Y_grid, z = t(lower),
        zlim = zlim_beta, col = cols_beta,
        xlab = paste("Grid points", predictor),
        ylab = paste("Grid points", x$info$call$formula[[2]]),
        main = paste0("Lower ", x$inference$conf_level * 100, "% Band")
      )

      # Plot upper confidence band
      image(
        x = X_grid, y = Y_grid, z = t(upper),
        zlim = zlim_beta, col = cols_beta,
        xlab = paste("Grid points", predictor),
        ylab = paste("Grid points", x$info$call$formula[[2]]),
        main = paste0("Upper ", x$inference$conf_level * 100, "% Band")
      )

      # Add color scale
      par(mar = legend_info$mar)
      fields::image.plot(
        legend.only = TRUE,
        zlim = zlim_beta,
        col = cols_beta,
        legend.mar = 4,
        legend.args = list(text = "Value", cex = 0.8)
      )
    } else {
      # Single panel plot for coefficient only
      fields::image.plot(
        x = X_grid, y = Y_grid, z = t(beta),
        zlim = zlim_beta, col = cols_beta,
        xlab = paste("Grid points", predictor),
        ylab = paste("Grid points", x$info$call$formula[[2]]),
        main = expression(paste("Estimated ", beta(r, s)))
      )
    }
  } else if (which == "t") {
    z <- x$inference$t_values[[predictor]][oY, oX, drop = FALSE]
    title <- expression(t - values)
    zlim <- max(abs(z), na.rm = TRUE) * c(-1, 1)
    cols <- colorRampPalette(c("blue", "white", "red"))(100)

    fields::image.plot(
      x = X_grid, y = Y_grid, z = t(z),
      zlim = zlim, col = cols,
      xlab = paste("Grid points", predictor),
      ylab = paste("Grid points", x$info$call$formula[[2]]),
      main = title
    )
  } else if (which == "p") {
    t_vals <- x$inference$t_values[[predictor]][oY, oX, drop = FALSE]
    z <- 2 * (1 - pt(abs(t_vals), df = x$stats$effective_df))
    title <- expression(p - values)

    # Set up plot parameters with extra margin on top for the legend
    par(mar = c(5, 4, 6, 4)) # Increased top margin

    # Create color scheme for p-values
    breaks <- c(0, 0.01, 0.05, 0.1, 1)
    cols <- colorRampPalette(c(
      "#800000", # dark red for highly significant (< 0.01)
      "#FF4500", # orangered for significant (0.01-0.05)
      "#FFA500", # orange for marginally significant (0.05-0.1)
      "#FFFACD" # light yellow for non-significant (> 0.1)
    ))(100)

    # Plot main image without legend
    image(
      x = X_grid, y = Y_grid, z = t(z),
      col = cols, zlim = c(0, 1), main = "p-values",
      xlab = paste("Grid points", predictor),
      ylab = paste("Grid points", x$info$call$formula[[2]])
    )

    # Add contour lines
    contour(
      x = X_grid, y = Y_grid, z = t(z),
      levels = c(0.01, 0.05, 0.1),
      add = TRUE, lty = c(1, 2, 3),
      lwd = 2, col = "darkblue",
      drawlabels = FALSE
    )

    # Add legend with fixed top-left positioning and smaller size
    par(xpd = TRUE) # Allow plotting outside the plot region
    usr <- par("usr") # Get the plot coordinates
    legend(
      x = usr[1], # Start at the left edge
      y = usr[4] + (usr[4] - usr[3]) * 0.15, # Position above the plot
      legend = c("p = 0.01", "p = 0.05", "p = 0.10"),
      lty = c(1, 2, 3),
      lwd = 2,
      col = "darkblue",
      horiz = TRUE,
      cex = 0.6, # Reduced from 0.8 to 0.6 for smaller text
      bty = "n"
    ) # No box around legend
    par(xpd = FALSE) # Reset to default
  }
}
