
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


train_ncl <- function(y_data, x_data_selected, lambda, incl_const = 1, 
                      l_rate = 0.5, n_iter = 100){
  
  
  num_models <- length(x_data_selected)
  num_vars <- ncol(x_data_selected[[1]])
  n_data <- length(y_data)
  
  # preset params and storage objects
  coef_matrix <- matrix(data = 0, nrow = num_models, ncol = num_vars + 1)
  gradient_coefs <- matrix(nrow = num_models, ncol = num_vars)
  prediction <- matrix(nrow = num_models, ncol = n_data)
  fc_prediction <- rep(0, times = n_data)
  
  # gd iteration
  for (i in c(1:n_iter)){
    
    # calculate predictions for all models:
    for (model in c(1:num_models)){
      prediction[model,] <- coef_matrix[model,1]*incl_const + 
        x_data_selected[[model]] %*% coef_matrix[model,-1]
      a <- x_data_selected[[model]] %*% coef_matrix[model,-1]
    }
    
    # calculate fc_prediction:
    fc_prediction <- colSums(prediction)/num_models
    
    
    # calculate gradients:
    gradient_const <- (rowSums(t(t(prediction) - y_data)) - 
                         lambda*( rowSums(t(t(prediction) - fc_prediction)) ))*incl_const
    
    for (model in c(1:num_models)){
      temp <- (prediction[model,] - y_data) - 
        lambda*(prediction[model,] - fc_prediction)
      temp <- temp %*% x_data_selected[[model]]
      gradient_coefs[model,] <- temp
    }
    
    # calc mean of gradient
    mean_grad_const <- gradient_const/n_data
    mean_grad_coefs <- gradient_coefs/n_data
    
    # recalculate coefs
    coef_matrix[,1] <- coef_matrix[,1] - l_rate*mean_grad_const
    coef_matrix[,-1] <- coef_matrix[,-1] - l_rate*mean_grad_coefs
    
    # recalculate fc_prediction
  }
  
  # return coef
  return(coef_matrix)
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
        coefs <- replace(coefs, is.na(coefs), 0)
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


forecast_ncl <- function(y_data, x_data_selected, init_train, incl_cons, lambda = 0.1){
  
  ### function that forecasts values in y_data with data from x_data_selected
  # uses expanding window, with initial training data length of "init_train"
  # estimated with ncl error, minimised by gradient descent.
  
  # inputs:
  # y_data: targets to forecast
  # x_data_selected: list of ind variables for each model
  # init_train: initial length of training data
  # incl_cons: == 1 if a constant is to be included, 0 otherwise
  # lambda: ncl penalty term
  
  # check if input has the expected properties, stop with error if not
  stopifnot((incl_cons == 0 | incl_cons == 1), class(init_train) == "numeric", 
            class(lambda) == "numeric")
  
  
  # get lengths/nums
  len <- length(y_data)
  num_models <- length(x_data_selected)
  num_vars <- ncol(x_data_selected[[1]])
  
  # preset return values and storages
  forecasted_vals <- matrix(nrow = num_models, ncol = len - init_train)
  coefs_estimate <- array(data = NA, dim = c(len - init_train, num_vars + 1, num_models))
  return_list <- list()
  resized_x_data_selected <- list()
  
  i <- 1
  
  for (window_length in c(init_train:(len-1))){
    
    # resite x_data_selected
    for (n in c(1:num_models)){
      resized_x_data_selected[[n]] <- x_data_selected[[n]][1:window_length,]
    }
    
    # train model on train data
    coefs <- train_ncl(y_data[1:window_length], resized_x_data_selected, 
                       lambda = lambda, incl_const = incl_cons)
    
    # get individual forecasts and coefs from trained models
    for (model in c(1:num_models)){
      forecasted_vals[model,i] <- predict_custom(x_data = x_data_selected[[model]][window_length+1,], 
                                                 coefs = coefs[model,])
      coefs_estimate[i,,model] <- coefs[model,]
    }
    
    i <- i + 1
  }
  
  
  # save ind forecasts and coefs in desired output format
  for (model in c(1:num_models)){
    temp_list <- list(forecasted_vals[model,], coefs_estimate[,,model])
    return_list[[model]] <- temp_list
  }
  
  return(return_list)
  
}


x_data_selector <- function(x_data, selector_matrix){
  
  ### select independent variables (columns) from x_data: 
  # "non-selected" columns are set to zero.
  
  # inputs: x_data is the matrix of all independent variables, size nxp
  # and selector_matrix is a matrix with 1s in the main diagonal for
  # columns to keep, size pxp
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(x_data) == "matrix", class(selector_matrix) == "matrix",
            ncol(x_data) == nrow(selector_matrix), ncol(x_data) == ncol(selector_matrix))
  
  # select desired data with matrix multiplication
  x_data_selected <- x_data %*% selector_matrix
  return(x_data_selected)
}


forecast_fc <- function(y_data, x_data, selector_mat_list, init_train, incl_cons, 
                        method = "traditional", lambda = 0.1){
  
  # forecasts y_data with a combination of forecasts using information in x_data.
  # individual models are determined by selector_mat_list
  
  # inputs:
  # y_data: data to be forecasted.
  # x_data: independent variables
  # selector_mat_list: list of n diagonal matrices. 
  # determines the individual models.
  # init_train: initial training window lenght.
  # incl_cons: if == 1, a constant is included.
  # method: determines if traditional or ncl training should be used.
  
  num_models <- length(selector_mat_list)
  num_ind_variables <- ncol(x_data)
  
  # preset variables
  x_data_selected <- list()
  ind_forecasts <- list()
  
  # list models
  for (i in c(1:num_models)){
    x_data_selected[[i]] <- x_data_selector(x_data, selector_mat_list[[i]])
  }
  
  # train forecasts
  
  if (method == "traditional"){
    for (i in c(1:num_models)){
      ind_forecasts[[i]] <- forecast(y_data, x_data_selected[[i]], init_train, 
                                     incl_cons, method = "analytic")
    }
  } else if (method == "ncl"){
    
    ind_forecasts <- forecast_ncl(y_data, x_data_selected, init_train, incl_cons, lambda = lambda)
    
  }
  
  
  ##  average forecasts
  # store individual forecasts and coefs in vars
  forecast_length <- length(ind_forecasts[[1]][[1]])
  
  forecasted_vals_ind <- matrix(nrow = num_models, ncol = forecast_length)
  coefs_ind <- array(dim = c(forecast_length, num_ind_variables + 1, num_models))
  
  for (i in c(1:num_models)){
    forecasted_vals_ind[i,] <- ind_forecasts[[i]][[1]]
    coefs_ind[,,i] <- ind_forecasts[[i]][[2]]
  }
  
  fc_forecasts <- apply(forecasted_vals_ind, 2, mean)
  fc_coefs <- apply(coefs_ind,c(1,2),mean)
  
  # return
  result_list <- list(fc_forecast = fc_forecasts, fc_coef = fc_coefs, 
                      ind_forecast = forecasted_vals_ind, ind_coef = coefs_ind)
  
  return(result_list)
  
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
X <- def_dgp_params(X, coefs = c(1,1,0,0), r_squared = 0.02, correl = 0.8, num_of_vars = 4, len = 200)
X <- simulate_dgp(X,1)

len <- 200
init_train <- 70

y_test_data <- X$DGP_1$sim_data$y[1:len]
x_test_data <- X$DGP_1$sim_data$x[1:len, 1:4]

forecast_test <- forecast(y_test_data, 
                          x_test_data, 
                          init_train, incl_cons = 0)

forecast_test_analytic <- forecast(y_test_data, x_test_data, init_train, 
                             incl_cons = 0, method = "analytic")

cw_test <- clark_west_test(y_test_data[(init_train+1):len], forecast_test[[1]], 
                forecast_hist_avg(y_test_data, init_train))

selector_mat_list <- list()

for (i in c(1:4)){
  
  diag_vect <- rep(0, times = 4)
  diag_vect[i] <- 1
  selector_mat_list[[i]] <- diag(diag_vect)  
}

trad_univariate_test <- forecast_fc(y_test_data, x_test_data, 
                                    selector_mat_list, init_train, incl_cons = 0)


### choose each univariate model

selector_mat_list <- list()

for (i in c(1:4)){
  
  diag_vect <- rep(0, times = 4)
  diag_vect[i] <- 1
  selector_mat_list[[i]] <- diag(diag_vect)  
}

x_data_selected <- list()

for (i in c(1:4)){
  x_data_selected[[i]] <- x_data_selector(x_test_data, selector_mat_list[[i]])
}

### test train_ncl

ncl_univariate_train_test <- train_ncl(y_test_data, x_data_selected, 0.1, incl_const = 0)


### test forecast_ncl

forecast_ncl(y_test_data, x_data_selected, init_train, 0)

ncl_univariate_test <- forecast_fc(y_test_data, x_test_data, 
                                    selector_mat_list, init_train, incl_cons = 0, method = "ncl", lambda = 1.0)

# get MSE
calc_MSE(trad_univariate_test$fc_forecast,y_test_data[71:200])
calc_MSE(ncl_univariate_test$fc_forecast,y_test_data[71:200])
calc_MSE(forecast_test[[1]],y_test_data[71:200])











### test

X <- list()
X <- def_dgp_params(X, coefs = c(1,1,0,0), r_squared = 0.05, correl = 0.6, num_of_vars = 4, len = 200)
X <- simulate_dgp(X,1)

len <- 200
init_train <- 70

y_test_data <- X$DGP_1$sim_data$y[1:len]
x_test_data <- X$DGP_1$sim_data$x[1:len, 1:4]

forecast_test <- forecast(y_test_data, 
                          x_test_data, 
                          init_train, incl_cons = 0)




selector_mat_list <- list()

for (i in c(1:4)){
  
  diag_vect <- rep(0, times = 4)
  diag_vect[i] <- 1
  selector_mat_list[[i]] <- diag(diag_vect)  
}

x_data_selected <- list()

for (i in c(1:4)){
  x_data_selected[[i]] <- x_data_selector(x_test_data, selector_mat_list[[i]])
}

trad_univariate_test <- forecast_fc(y_test_data, x_test_data, 
                                    selector_mat_list, init_train, incl_cons = 0)

ncl_univariate_test <- forecast_fc(y_test_data, x_test_data, 
                                   selector_mat_list, init_train, incl_cons = 0, method = "ncl", lambda = 0.2)

calc_MSE(trad_univariate_test$fc_forecast,y_test_data[71:200])
calc_MSE(ncl_univariate_test$fc_forecast,y_test_data[71:200])
calc_MSE(forecast_test[[1]],y_test_data[71:200])


