test_that("forecast_u_CAViaR runs without errors and returns expected structure", {
  # Generate sample data
  set.seed(123)
  #dates <- seq.Date(from = as.Date("2023-01-01"), by = "day", length.out = 1000)
  #asset_df <- data.frame(Date = dates, Return = rnorm(1000))
  df <- read.csv("C:/Users/chris/RStudioProjects/caviar/data/clean_returns.csv", row.names = 1)
  asset_df <- df["DEPYc1"]
  asset_df$Date <- as.Date(rownames(asset_df))
  colnames(asset_df) <- c("Return", "Date")
  rownames(asset_df) <- NULL

  print(tail(asset_df, 10))

  options(warn = 1)

  # Parameters
  c <- 0.05
  n <- 989
  m <- 250
  r <- 10
  verbose <- 3
  control.optim <- list(factr = 1e-7)
  control.solnp <- list(tol = 1e-7, trace = 0)
  var_model <- "AS"
  es_model <- "MULT"

  # Call the function
  result <- RollCAViaR(asset_df, c = c, n = n, m = m, r = r,
                       var_model = var_model, es_model = es_model,
                       verbose = verbose, control.optim = control.optim,
                       control.solnp = control.solnp)

  print(result)

  # Check the result is of the correct type
  expect_s3_class(result, "xts")

  # Check the result has the correct dimensions
  expect_equal(nrow(result), n)  # Ensure it has 'n' rows
  expect_equal(ncol(result), 2 * length(c))  # Should have columns for VaR and ES

  # Check column names
  expect_named(result, c("VaR", "ES"))

  # Validate the data range
  #expect_true(all(result$VaR >= 0))  # VaR should be negative
  #expect_true(all(result$ES >= 0))   # ES should be negative
})
