
# ffr: Functional Factor Regression

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/ffr)](https://CRAN.R-project.org/package=ffr)
[![R-CMD-check](https://github.com/luiswn/ffr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luiswn/ffr/actions)
<!-- badges: end -->

The `ffr` package provides a comprehensive framework for estimation, prediction and statistical inference in functional linear models. It implements the functional factor regression approach for modeling relationships between functional and scalar predictors and functional responses, originally introduced by Otto & Winter (2025).

## Preprint

arxiv link

## Features

- **Function-on-Function Regression**: Fit linear models with functional responses and predictors
- **Factor Estimation**: Automatically determine the correct number of factors using the eigenvalue difference method
- **Statistical Inference**: Generate confidence regions for coefficient surfaces
- **Data Preprocessing**: Transform irregularly spaced functional data onto consistent grids
- **Visualization**: Plot coefficient surfaces, p values, and more

## Installation

You can install the development version of `ffr` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("luiswn/ffr")
```

## Example

This is a minimal example to demonstrate the core functionalities of the `ffr` package. The code example does not coincide with the application in the paper.

### Data Preparation

The data used in the following stems from the [epftoolbox](https://epftoolbox.readthedocs.io/en/latest/index.html). Please install the Python package as explained [here](https://epftoolbox.readthedocs.io/en/latest/modules/started.html#installation) before you execute the R code.

``` r
library(reticulate)
library(tidyverse)

# In case of multiple Python installations, set path to use_python("/path/to/python/epftoolbox")

# Import python "epftoolbox" and read the data
py_run_string("
import epftoolbox.data
result = epftoolbox.data.read_data(path='.', dataset='DE', years_test=0)
df = result[0]                # Extract first element from tuple
df_unique = df.reset_index()  # Reset index to handle duplicate timestamps
")

# Reshape electricity data and prepare variables
DE_epf_data <- as.data.frame(py$df_unique) %>%
  mutate(index = as.character(index)) %>%
  separate(index, c("Day", "Time"), "\\s+") %>%
  rename("Load" = `Exogenous 1`,
         "Wind_Solar_Gen" = `Exogenous 2`) %>%
  mutate(Time = rep(0:23, times=2184),
         Day = rep(1:2184,each=24),
         Load = Load / 1000,
         Wind_Solar_Gen =  Wind_Solar_Gen / 1000)

reshape_to_wide <- function(data, value_column) {
  pivot_wider(
    data %>% select(Day, Time, !!sym(value_column)),
    names_from = Time,
    values_from = !!sym(value_column)
  ) %>%
    select(-Day)
}

Price <- reshape_to_wide(DE_epf_data, "Price")
Load <- reshape_to_wide(DE_epf_data, "Load")
Wind_Solar_Gen <- reshape_to_wide(DE_epf_data, "Wind_Solar_Gen")

# Concatenate regression variables to data list
data <- list(Price = as.matrix(Price), 
             Load = as.matrix(Load), 
             Wind_Solar_Gen = as.matrix(Wind_Solar_Gen))
```
### Factor Estimation

``` r
library(ffr)

# Estimate the number of factors for all functional predictor
K_EDresult <- fed(Price ~ . , data, gamma = 93)
print(K_EDresult$K)

```

### Fit Fully Functional Linear Regression

``` r
# Fit a function-on-function linear model with scalar predictor
ffr_model <- flm(Price ~ . , data, K = K_EDresult$K)

# Save model summary
ffr_model_sum <- summary(ffr_model)

# Print model summary
print(ffr_model_sum)

```

### Visualization

``` r
# Visualize the coefficient and confidence region heatmaps for the electricity load predictor
plot(ffr_model_sum, predictor = "Load", which = "beta", conf.region = TRUE)

# Visualize the 3D coeffcient plot for the solar and windgeneration predictor
plot(ffr_model_sum, predictor = "Wind_Solar_Gen", which = "beta_3D", conf.region = FALSE)

# Visualize the pointwise p value heatmap for the electricity load predictor
plot(ffr_model_sum, predictor = "Load", which = "p")

# View functional R² across the response domain
plot(ffr_model_sum, which = "R2")

```

## Citation

If you use this package in a scientific context, please cite:

Otto, S., & Winter, L. (2025). Functional factor regression with an application to electricity price curve modeling.

## License

MIT
