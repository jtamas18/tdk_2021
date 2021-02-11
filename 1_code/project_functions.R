
def_dgp_params <- function(X, coefs, r_squared, correl, num_of_vars = 14,
                               len = 300, sim_runs = 100, dgp_name = ""){
  
  # stores the data generating process parameters in list "dgp_param".
  # returns list "X" with dgp_param" as first elemenet.
  # X will store all other values related to the analysis
  
  # inputs: 
  # storage list X
  # DGP parameters
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(X) == "list", class(coefs) == "numeric", class(r_squared) == "numeric",
            class(correl) == "numeric", class(num_of_vars) == "numeric", 
            class(len) == "numeric", class(sim_runs) == "numeric", 
            class(dgp_name) == "character", length(coefs) == num_of_vars,
            r_squared >= 0, r_squared <= 1, correl >= -1, correl <= 1)
  
  # store parameters:
  dgp_case <- list()      
  dgp_param <- list()
  dgp_param$coefs <- coefs
  dgp_param$r_squared <- r_squared
  dgp_param$correl <- correl
  dgp_param$num_of_vars <- num_of_vars
  dgp_param$ts_length <- len
  dgp_param$sim_runs <- sim_runs
  
  dgp_param$resized_coefs <- ""
  dgp_param$covar_matrix <- ""
  
  dgp_case$dgp_param <- dgp_param
  
  # append DGP list to storage matrix X:
  X[[length(X)+1]] <- dgp_case
  
  # call the num DGP "DGG_dgp_name" if dgp_name is specified, call "DGP_[n]" otherwise
  if (dgp_name == ""){
    names(X) <- c(names(X[1:length(names(X))-1]), paste0("DGP_", toString(length(X))))
  } else {
    names(X) <- c(names(X[1:length(names(X))-1]), paste0("DGP_", dgp_name))
  } 
  
  return(X)
}




simulate_dgp <- function(X, dgp_index){
  # simulates the DGP defined in "dgp_parameters" times "dgp_parameters$sim_num"
  # stores the sim_num different trajectories in a matrix (size= ts_length * sim_num)
  # where each column is a different trajecotry
  # appends the matrix to X$DGP_[n]$trajectories  
  
  # inputs: 
  # the storage list X
  # the index of the DGP in X; eg. for dgp_index = 1, values used in the simulation
  # are taken from "X$DGP_1$..."
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(X) == "list", class(dgp_index) == "numeric")
  
  
  # take parameters from "X$DGP_[n]$dgp_parameters:
  ts_length <- X[[dgp_index]]$dgp_param$ts_length
  coefs <- X[[dgp_index]]$dgp_param$coefs
  r_squared <- X[[dgp_index]]$dgp_param$r_squared
  correl <- X[[dgp_index]]$dgp_param$correl
  num_of_vars <- X[[dgp_index]]$dgp_param$num_of_vars
  sim_runs <- X[[dgp_index]]$dgp_param$sim_runs
  
  ### generate the x_i-s:
  
  covar_first_row <- rep(1, times = num_of_vars) 
  for (i in c(2:num_of_vars)){
    covar_first_row[i] <- correl^(i-1)
  }
  
  covar_matrix <- toeplitz(covar_first_row)
  
  
  x <- mvtnorm::rmvnorm(ts_length*sim_runs, mean = rep(0, times = num_of_vars), 
                        sigma = covar_matrix )
  
  ### resize the coefs to match the desired r_squared:
  
  if (sum(coefs) == num_of_vars){
    # equal predictive power case
    denom <- r_squared
    nomin <- (1 - r_squared)*sum(covar_matrix)
    resized_coefs <- sqrt(denom/nomin)*coefs
  } else {
    # "sparse" case
    denom <- r_squared
    
    covar_matrix_temp <- covar_matrix
    covar_matrix_temp[row(covar_matrix_temp) > num_of_vars/2 | 
                        col(covar_matrix_temp) > num_of_vars/2] <- 0
    
    nomin <- (1 - r_squared)*sum(covar_matrix_temp)
    resized_coefs <- sqrt(denom/nomin)*coefs
  }
  
  
  ### generate y_i-s:
  
  y <- rowSums(t(resized_coefs * t(x))) + rnorm(ts_length*sim_runs, mean = 0, sd = 1)
  
  # store generated time series in storage matrix X
  X[[dgp_index]]$sim_data$x <- x
  X[[dgp_index]]$sim_data$y <- y
  
  # store new params in dgp_param
  X[[dgp_index]]$dgp_param$resized_coefs <- resized_coefs
  X[[dgp_index]]$dgp_param$covar_matrix <- covar_matrix
  
  return(X)
  
}


train_gradient <- function(y_data, x_data, incl_const = 1, 
                              l_rate = 0.5, n_iter = 100){
  
  ### traines a linear model with gradient descent. 
  # returns the estiamted coefficients in a vector.
  # first element of the vector is the constant (even if set to zero),
  # other elements are the slopes of the independent variables in order.
  
  # inputs:
  # y_data: targets used for training. numeric vector/matrix.
  # x_data: numeric matrix of predictors. nrow must be equal to length og y_data.
  # incl_const: numeric scalar, either 1 or 0. If == 1, a constant is estimated, 
  # if == 0, constant is set to 0.
  # l_rate is the learning rate, 
  # n_iter is the number of iterations in the gradient descent.
  
  
  
  # check if input has the expected properties, stop with error if not
  stopifnot((incl_const == 0 | incl_const == 1), class(l_rate) == "numeric", 
            class(n_iter) == "numeric", nrow(x_data) == length(y_data))
  
  # preset params
  const <- 0
  b_hat <- rep(0, times = ncol(x_data))
  
  x_data <- t(x_data)
  n_data <- length(y_data)
  
  for (i in c(1:n_iter)){
    # calculate prediction
    prediction <- const*incl_const + t(b_hat) %*% x_data
    
    # calculate gradients
    
    gradient_const <- sum((prediction - y_data)*incl_const)
    gradient_b_hat <- x_data %*% t(prediction - y_data)
    
    # calculate mean of gradient
    mean_grad_const <- gradient_const/n_data 
    mean_grad_b_hat <- gradient_b_hat/n_data 
    
    # update parameters by small part of averaged gradient
    const <- const - l_rate * mean_grad_const 
    b_hat <- b_hat - l_rate * mean_grad_b_hat
  }
  
  return(c(const, b_hat))
}

predict_custom <- function(x_data, coefs){
  
  ### makes a prediction based on the linear model with coefs
  # and independent variables in x_data. returns prediction
  
  prediction <- coefs[1] + x_data %*% coefs[-1]
  
  return(prediction)
} 

forecast <- function(y_data,x_data,init_train, incl_cons, method = "gd"){
  
  ### function that forecasts values in y_data with variables from x_data.
  # uses expanding window, with initial training data length of "init_train"
  # can be estimated by "analytic" or "gradient descent" (gd) methods.
  
  # inputs:
  # y_data: targets to forecast
  # x_data: predictors
  # init_train: initial length of training data
  # incl_cons: == 1 if a constant is to be included, 0 otherwise
  # method: == "gd" if estimation is by gradient descent, "analytic" 
  # if estimation by traditinal matrix inversion.
  
  # check if input has the expected properties, stop with error if not
  stopifnot((incl_cons == 0 | incl_cons == 1), class(init_train) == "numeric", 
            (method == "gd" | method == "analytic"), nrow(x_data) == length(y_data))
  
  # get length of data
  len <- length(y_data)
  
  # preset return values
  forecasted_vals <- rep(NA, times = len - init_train)
  coefs_estimate <- matrix(data = NA, nrow = len - init_train, 
                           ncol = ncol(x_data) + 1)
  i <- 1
  
  if (method == "gd"){
    for (window_length in c(init_train:(len-1))){
    
      # train model on train data
      coefs <- train_gradient(y_data[1:window_length], x_data[1:window_length,], 
                   incl_const = incl_cons)
    
      # forecast with trained model
      pred <- predict_custom(x_data = x_data[window_length+1,], coefs = coefs)
    
      # save forecast and estimated coefs
      forecasted_vals[i] <- pred
      coefs_estimate[i,] <- coefs
    
      i <- i + 1
    }
  } else if (method == "analytic"){
    
    for (window_length in c(init_train:(len-1))){
      
      # train lin model by matrix inversion
      # do not include intercept if incl_cons == 0
      if (incl_cons == 1){
        lin_model <- lm(formula = y_data[1:window_length] ~ ., 
                        data = as.data.frame(x_data[1:window_length,]))
        # get coefs
        coefs <- lin_model$coefficients
      } else {
        lin_model <- lm(formula = y_data[1:window_length] ~ . + 0, 
                        data = as.data.frame(x_data[1:window_length,]))
        # get coefs
        coefs <- c(0,as.vector(lin_model$coefficients))
      }
      
      
      # forecast with trained model
      pred <- predict_custom(x_data = x_data[window_length + 1,], coefs = coefs)
      
      # save forecast and estimated coefs
      forecasted_vals[i] <- pred
      coefs_estimate[i,] <- coefs
      
      i <- i + 1
      }
  }
  
  
  return(list(forecasted_vals, coefs_estimate))
    
}

forecast_hist_avg <- function(y_data,init_train){
  
  ### calculates forecasts in y_data based on previous historical average.
  ### uses expnading window with initial training window length of "init_train".
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(y_data) == "numeric", class(init_train) == "numeric", 
            length(y_data) > init_train)
  
  
  # get length of data
  len <- length(y_data)
  
  # preset return value
  forecasted_vals <- rep(NA, times = len - init_train)
  i <- 1
  
  
  for (window_length in c(init_train:(len-1))){
    
    # train lin model with only intercept
    hist_avg <- mean(y_data[1:window_length])
    
    # save hist_avg as forecast
    forecasted_vals[i] <- hist_avg
    
    i <- i + 1
  }
  
  return(forecasted_vals)
}

  
calc_MSE <- function(predictions, targets){
  
  ### calculates MSE, with "predictons" as a vector of predictions 
  ### and "targets" as a vector of targets.
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(predictions) == "numeric", class(targets) == "numeric",
            length(predictions) == length(targets))
  
  # calculate MSE
  MSE <- sum((predictions - targets)^2)*(1/length(predictions))
  
  return(MSE)
}

calc_r_squared <- function(targets, predictions, hist_avg_pred){
  
  ### calculates out of sample R^2 of the predictions in vector "predictions".
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(predictions) == "numeric", class(targets) == "numeric", 
            class(hist_avg_pred) == "numeric", length(predictions) == length(targets), 
            length(predictions) == length(targets))
  
  MSE_pred <- calc_MSE(predictions, targets)
  MSE_hist_avg <- calc_MSE(hist_avg_pred, targets)
  r_squared <- 1 - MSE_pred/MSE_hist_avg
  
  return(r_squared)
}

clark_west_test <- function(targets, predictions, hist_avg_pred){
  
  ### calculates t-statistic and p-values for one-sided Clark & West test
  # null hypothesis: no difference between predictive power of predictions 
  # and benchmark. Alt hypothesis: predictions better forecasts targets 
  # than benchmark.
  
  # calculates t-  & p-values as in Rapach(2009).

  # check if input has the expected properties, stop with error if not
  stopifnot(class(predictions) == "numeric", class(targets) == "numeric",
            class(hist_avg_pred) == "numeric", length(predictions) == length(targets),
            length(hist_avg_pred) == length(targets))
  
  # calculate data for clark & west regression
  c_w_data <- (targets - hist_avg_pred)^2 - ((targets - predictions)^2 - (hist_avg_pred - predictions)^2)
  
  c_w_model <- lm(formula = c_w_data ~ 1)
  
  t_stat <- summary(c_w_model)[["coefficients"]][, "t value"]
  p_value <- 1 - pnorm(t_stat)
  
  return(c(t_stat, p_value))
  
}
  
### próba
X <- list()
X <- def_dgp_params(X, coefs = c(1,1,0,0), r_squared = 0.9, correl = 0.3, num_of_vars = 4, len = 500)
X <- simulate_dgp(X,1)

len <- 200
init_train <- 70

y_test_data <- X$DGP_1$sim_data$y[1:len]
x_test_data <- X$DGP_1$sim_data$x[1:len, 1:4]

forecast_test <- forecast(y_test_data, 
                          x_test_data, 
                          init_train, incl_cons = 0)
forecast_test_an <- forecast(y_test_data, x_test_data, init_train, 
                             incl_cons = 0, method = "analytic")

cw_test <- clark_west_test(y_test_data[(init_train+1):len], forecast_test[[1]], 
                forecast_hist_avg(y_test_data, init_train))
cw_test
