#' Rolling window forecast of VaR and ES using CAViaR
#'
#'
#' @useDynLib caviar
#' @export
#' @import progress
#' @importFrom xts xts
#' @importFrom zoo na.locf
#' @importFrom solnp solnp
RollCAViaR <- function(
    data,
    c,
    n,
    m,
    r = 10,
    var_model = "ADAPTIVE",
    es_model = "MULT",
    loss = "AL",
    lb = -1,
    ub = 1,
    grid_size = 10000,
    verbose = 0,
    ...
) {



  # Initialize variables for es_model
  if (es_model == "MULT") {
    es_func <- mult_ES_C
    n_gammas <- 1
    upper_gamma <- ub
    lower_gamma <- lb
  } else if (es_model == "AR") {
    es_func <- ar_ES_C
    n_gammas <- 3
    upper_gamma <- ub
    lower_gamma <- 0
  } else {
    stop("Invalid ES model")
  }

  # Initialize variables for var_model
  if (var_model == "ADAPTIVE") {
    var_func <- adaptive_C
    n_betas <- 1
    lower_beta <- 0
  } else if (var_model == "SAV") {
    var_func <- SAV_C
    n_betas <- 3
    upper_beta <- c(ub, 1, ub)
    lower_beta <- c(0, 0, 0)
  } else if (var_model == "AS") {
    var_func <- AS_C
    n_betas <- 4
    upper_beta <- c(1, ub, ub, ub)
    lower_beta <- c(0, 0, 0, 0)
  } else if (var_model == "indirectGARCH") {
    var_func <- indirectGARCH_C
    n_betas <- 3
    upper_beta <- c(ub, ub, ub)
    lower_beta <- c(0, 0, 0)
  } else if (var_model == "linearGARCH") {
    var_func <- linearGARCH_C
    n_betas <- 3
    upper_beta <- c(ub, ub, ub)
    lower_beta <- c(lb, lb, lb)
  } else if (var_model == "linearTGARCH") {
    var_func <- linearTGARCH_C
    n_betas <- 4
    upper_beta <- c(ub, ub, ub, ub)
    lower_beta <- c(lb, lb, lb, lb)
  } else if (var_model == "gjrGARCH") {
    var_func <- gjrGARCH_C
    n_betas <- 4
    upper_beta <- c(ub, ub, ub, ub)
    lower_beta <- c(lb, lb, lb, lb)
  } else {
    stop("Invalid VaR model")
  }

  if (loss == "AL") {
    loss_func <- al_log_loss_function_C
  } else {
    stop("Invalid loss function")
  }

  # Record total number of parameters
  n_params <- n_betas + n_gammas

  # Merge constraints
  lower <- c(lower_beta, lower_gamma)
  upper <- c(upper_beta, upper_gamma)

  # defining progress bar
  if (verbose > 0) {
    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsed || Estimated time remaining: :eta]",
      total = n,
      complete = "=",
      incomplete = "-",
      current = ">",
      show_after = 0.2,
      clear = TRUE
      )
  }

  # Initialize variables before iterating for rolling window
  df <- tail(data$Return, n + m)
  var <- matrix(NA, nrow = n, ncol = length(c))
  es <- matrix(NA, nrow = n, ncol = length(c))
  last_params <- NA
  for (i in 1:n) {
    if (verbose > 0) pb$tick()
    # Retrieve the rolling window
    window_start <- i
    window_end <- m + i - 1
    window <- df[window_start:window_end]

    # Reparametrize for every r period
    if (i %% r == 0 || i == 1) {
      optim_function <- objective_handler(window, c, var_func, es_func, n_betas, loss_func)

      # Initialize grid for parameter search
      grid <- matrix(runif(grid_size * n_params, min = lower, max = upper), ncol = n_params)
      grid_losses <- apply(grid, 1, optim_function)
      best_grid_index <- which.min(grid_losses)
      initial_params <- grid[best_grid_index, ]

      # Optimize using optim with initial parameters from grid
      result <- tryCatch({
        solnp(
          par = initial_params,
          fn = optim_function,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          ...
        )
      }, error = function(e) {
        message("Optimization failed: ", e$message)
        NULL
      }, warning = function(w) {
        message("Optimization warning: ", w$message)
        NULL
      })

      # Update parameters
      if (!is.null(result$par)) {
        last_params <- result$par

        if (result$convergence != 0) {
          warning("Optimization did not converrge iteration ", i, ", convergence code: ", result$convergence, ", message: ", result$message)
        }
      } else {
        warning("Optimization failed at iteration ", i, ", skipping param optimization")
        if(any(is.na(last_params))) {
          var[i] <- NA
          es[i] <- NA
          next
        }
      }
    }

    # Get parameters
    betas <- last_params[1:n_betas]
    gammas <- last_params[n_betas + 1:n_params]

    # Forecast VaR for next timestep
    u <- quantile(window, probs = c)
    var[i] <- var_func(y = window, betas = betas, u = u, c = c)[m]

    # Initialize vectors for current and next timestep for y, x, and Q - which is the empirical quantile first iteration
    Q_input <- if (i == 1) c(u, var[i]) else var[(i - 1):i]
    y_input <- c(window[m], 0)
    if (es_model == "ar") {
      x_input <- if (i == 1) rep(0, 2) else c(var[i - 1] - es[i - 1], 0)
      es_output <- es_func(y = y_input, Q = Q_input, gammas = gammas, x = x_input)
    } else {
      es_output <- es_func(y = y_input, Q = Q_input, gammas = gammas)
    }
    es[i] <- es_output[2]
  }

  if (verbose > 0) {
    pb$terminate(); invisible()
  }

  # Create xts objects for VaR and ES
  dates <- tail(data$Date, n)
  VaR <- xts(-as.numeric(var), order.by = dates, colnames = "VaR")
  ES <- xts(-as.numeric(es), order.by = dates, colnames = "ES")
  results_xts <- merge(VaR, ES)

  if (any(is.na(results_xts))) {
    warning("There were NA values, carry forward last non-NA")
    results_xts <- zoo::na.locf(results_xts, na.rm = FALSE)
  }

  return(results_xts)
}

objective_handler <- function(y, c, func_Q, func_ES, n_betas, loss_function) {
  function(params) {
    betas <- params[1:n_betas]
    gammas <- params[(n_betas + 1):length(params)]

    u <- quantile(y, probs = c)

    Q <- func_Q(y = y, betas = betas, u = u, c = c)
    ES <- func_ES(y = y, Q = Q, gammas = gammas)

    # Calculate the loss vector
    ES[ES >= 0] <- -1e-10
    loss_vector <- loss_function(y, Q, ES, c)

    out_of_bounds_val <- 5 * max(loss_vector)

    # # Handle NA or infinite values to prevent errors in optimization
    loss_vector[is.na(loss_vector) | is.infinite(loss_vector)] <- out_of_bounds_val
    loss_vector[Q >= 0] <- out_of_bounds_val

    return(mean(loss_vector))
  }
}
