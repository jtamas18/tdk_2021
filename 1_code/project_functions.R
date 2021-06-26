
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




simulate_dgp <- function(X, dgp_index, covar = "toeplitz"){
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
  
  if (covar == "toeplitz"){
    covar_first_row <- rep(1, times = num_of_vars) 
    for (i in c(2:num_of_vars)){
      covar_first_row[i] <- correl^(i-1)
    }
    
    covar_matrix <- toeplitz(covar_first_row)
  } else if (covar == "constant"){
    covar_matrix <- matrix(data = correl, nrow = num_of_vars, ncol = num_of_vars)
    
    for (i in c(1:num_of_vars)){
      covar_matrix[i,i] <- 1
    }
  }
  
  
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
                              l_rate = 0.1, n_iter = 100){
  
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
                      l_rate = 0.5, n_iter = 25, OLS_coefs = NA){
  
  
  num_models <- length(x_data_selected)
  num_vars <- ncol(x_data_selected[[1]])
  n_data <- length(y_data)
  
  # preset params and storage objects
  coef_matrix <- matrix(data = 0, nrow = num_models, ncol = num_vars + 1)
  gradient_coefs <- matrix(nrow = num_models, ncol = num_vars)
  prediction <- matrix(nrow = num_models, ncol = n_data)
  fc_prediction <- rep(0, times = n_data)
  
  # preset starting coefs at OLS coefs if provided
  if (is.numeric(OLS_coefs)){
    coef_matrix <- t(OLS_coefs)
  }
  
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

forecast <- function(y_data,x_data,init_train, incl_cons, method = "gd", l_rate = 0.5, n_iter = 25){
  
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
                   incl_const = incl_cons, l_rate = l_rate, n_iter = n_iter)
    
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
        coefs <- replace(coefs, is.na(coefs), 0)
      } else {
        lin_model <- lm(formula = y_data[1:window_length] ~ . + 0, 
                        data = as.data.frame(x_data[1:window_length,]))
        # get coefs
        coefs <- c(0,as.vector(lin_model$coefficients))
        coefs <- replace(coefs, is.na(coefs), 0)
      }
      
      # does not do anything
      l_rate <- l_rate + 1
      n_iter <- n_iter + 1
      
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


forecast_ncl <- function(y_data, x_data_selected, init_train, incl_cons, lambda = 0.1, 
                         l_rate = 0.5, n_iter = 25, OLS_coefs = NA){
  
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
                       lambda = lambda, incl_const = incl_cons, l_rate = l_rate, 
                       n_iter = n_iter, OLS_coefs = OLS_coefs[i,,])
    
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



forecast_pncl <- function(y_data, x_data_selected, init_train, incl_cons, 
                          lambda = 0.1, l_rate = 0.5, n_iter = 25){
  
  ### function that forecasts values in y_data with data from x_data_selected
  # uses expanding window, with initial training data length of "init_train"
  # estimated with ncl error, minimised by gradient descent.
  
  # inputs:
  # y_data: targets to forecast
  # x_data_selected: list of ind variables for each model
  # init_train: initial length of training data
  # incl_cons: == 1 if a constant is to be included, 0 otherwise
  # lambda: ncl penalty term
  # noNCL index: index of model in data selected that should not be learn with NCL
  
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
    
    # preset coefs storage
    coefs <- matrix(data = NA, nrow = num_models, ncol = num_vars + 1)
    
    # train NCL models on train data
    coefs[1:(num_models-1),] <- train_ncl(y_data[1:window_length], resized_x_data_selected[-num_models], 
                                          lambda = lambda, incl_const = incl_cons, 
                                          l_rate = l_rate, n_iter = n_iter)
    
    
    
    # train lin model by matrix inversion
    # do not include intercept if incl_cons == 0
    if (incl_cons == 1){
      lin_model <- lm(formula = y_data[1:window_length] ~ ., 
                      data = as.data.frame(resized_x_data_selected[num_models]))
      # get coefs
      temp_coefs <- lin_model$coefficients
    } else {
      lin_model <- lm(formula = y_data[1:window_length] ~ . + 0, 
                      data = as.data.frame(resized_x_data_selected[num_models]))
      # get coefs
      temp_coefs <- c(0,as.vector(lin_model$coefficients))
    }
    temp_coefs <- replace(temp_coefs, is.na(temp_coefs), 0)
    
    
    
    # train trad model with OLS
    coefs[num_models,] <- temp_coefs
    
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
                        method = "traditional", lambda = 0.1, l_rate = 0.5, n_iter = 25, OLS_coefs = NA){
  
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
                                     incl_cons, method = "analytic", l_rate = l_rate, n_iter = n_iter)
      
    }
  } else if (method == "ncl"){
    
    ind_forecasts <- forecast_ncl(y_data, x_data_selected, init_train, incl_cons, 
                                  lambda = lambda, l_rate = l_rate, n_iter = n_iter, OLS_coefs = OLS_coefs)
    
  } else if (method == "p_ncl"){
    ind_forecasts <- forecast_pncl(y_data, x_data_selected, init_train, incl_cons, lambda = lambda,
                                   l_rate = l_rate, n_iter = n_iter)
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
  #stopifnot(class(predictions) == "numeric", class(targets) == "numeric",
  #          length(predictions) == length(targets))
  
  # calculate MSE
  MSE <- sum((predictions - targets)^2)*(1/length(predictions))
  
  return(MSE)
}

calc_r_squared <- function(targets, predictions, hist_avg_pred){
  
  ### calculates out of sample R^2 of the predictions in vector "predictions".
  
  # check if input has the expected properties, stop with error if not
  #stopifnot(class(targets) == "numeric", 
  #          class(hist_avg_pred) == "numeric", length(predictions) == length(targets))
  
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
  
  # calculates t  & p-values as in Rapach(2009).

  # check if input has the expected properties, stop with error if not
  stopifnot(class(predictions) == "numeric", class(targets) == "numeric",
            class(hist_avg_pred) == "numeric", length(hist_avg_pred) == length(targets))
  
  # calculate data for clark & west regression
  c_w_data <- (targets - hist_avg_pred)^2 - ((targets - predictions)^2 - (hist_avg_pred - predictions)^2)
  
  c_w_model <- lm(formula = c_w_data ~ 1)
  
  t_stat <- summary(c_w_model)[["coefficients"]][, "t value"]
  p_value <- 1 - pnorm(t_stat)
  
  return(c(t_stat, p_value))
  
}

choose_models <- function(x_data, method = "univariate", k = NA){
  
  # creates list of model selector matrixes
  # input: x_data: matrix of predictors, method == "univariate" or "csr" 
  # k is only needed if method == "csr"
  
  # check if input has the expected properties, stop with error if not
  stopifnot(class(x_data) == "matrix", class(method) == "character", (class(k) == "numeric" | is.na(k)))
  
  # preset return value
  selector_mat_list <- list()
  
  # get num of models
  num_vars <- ncol(x_data)
  
  if (method == "univariate"){
    
    num_models <- num_vars
    
    # create selector matrixes
    for (i in c(1:num_models)){
      
      diag_vect <- rep(0, times = num_vars)
      diag_vect[i] <- 1
      selector_mat_list[[i]] <- diag(diag_vect)  
    }
    
  } else if (method == "csr"){
   num_models <- choose(num_vars,k)
   
   # create selector matrixes
   #### NOT IMPLEMENTED
   for (i in c(1:num_models)){
     
     diag_vect <- rep(0, times = num_vars)
     diag_vect[i] <- 1
     selector_mat_list[[i]] <- diag(diag_vect)  
   } 
   print("not writen yet! does not work as intended!")
  } else if (method == "kitchen_sink"){
    
    #create selector matrix (identity matrix)
    diag_vect <- rep(1, times = num_vars)
    selector_mat_list[[1]] <-  diag(diag_vect)
    
  }
  
  # return list of selector matrixes
  return(selector_mat_list)
}


sim_forecasts <- function(X_dgps, init_train = 70, method = "traditional" , 
                          models = "univariate", k = NULL, lambda = NULL){
  
  tic()
  # check if input has the expected properties, stop with error if not
  stopifnot(class(X_dgps) == "list", class(init_train) == "numeric",
            class(models) == "character")
  
  # get num of dgps
  num_dgps <- length(X_dgps)
  
  # iterate over dgps
  for (dgp in c(1:num_dgps)){
    
    # get dgp params
    ts_length <- X_dgps[[dgp]]$dgp_param$ts_length
    sim_runs <- X_dgps[[dgp]]$dgp_param$sim_runs
    y_data <- X_dgps[[dgp]]$sim_data$y
    x_data <- X_dgps[[dgp]]$sim_data$x
    
    # preset data start and end indices
    data_start <- 1
    data_end <- ts_length
    
    # create variable selector matrices
    selector_mat_list <- choose_models(x_data = x_data, method = models, k = k)
    
    # iterate over sim nums
    for (sim_run in c(1:sim_runs)){
      
      # get data for this "run":
      y_this_run <- y_data[data_start:data_end]
      x_this_run <- x_data[data_start:data_end,]
      
      # calc forecasts
      forecast_results <- forecast_fc(y_data = y_this_run, x_data = x_this_run, 
                                      selector_mat_list, init_train, incl_cons = 0, method = method, lambda = lambda)
      
      # save forecasts baed on model & method type
      X_dgps[[dgp]]$forecast_results$temp_name[[sim_run]] <- forecast_results
      
      # go to next sim_run:
      data_start <- data_start + ts_length
      data_end <- data_end + ts_length
      
    }
    
    # rename forecast results list
    new_name <- paste0(method, "_", models)
    if (method == "ncl" | method == "P_ncl"){
      new_name <- paste0(new_name, "_lambda_", as.character(lambda))
    }
    if (models == "csr"){
      new_name <- paste0(new_name, "_k_", as.character(k))
    }
    
    forecast_num <- length(X_dgps[[dgp]]$forecast_results)
    names(X_dgps[[dgp]]$forecast_results) <- c(names(X_dgps[[dgp]]$forecast_results)[-forecast_num], 
                                               new_name)
    
  }
  toc()
  return(X_dgps)
}





sim_hist_avg_fc <- function(X_dgps, init_train){
  # get num of dgps
  num_dgps <- length(X_dgps)
  
  # iterate over dgps
  for (dgp in c(1:num_dgps)){
    
    # get dgp params
    ts_length <- X_dgps[[dgp]]$dgp_param$ts_length
    sim_runs <- X_dgps[[dgp]]$dgp_param$sim_runs
    y_data <- X_dgps[[dgp]]$sim_data$y
    
    # preset data start and end indices
    data_start <- 1
    data_end <- ts_length
    
    # iterate over sim nums
    for (sim_run in c(1:sim_runs)){
      
      # get data for this "run":
      y_this_run <- y_data[data_start:data_end]
      
      # calc forecasts
      forecast_results <- forecast_hist_avg(y_data = y_this_run, init_train)
      
      # save forecasts baed on model & method type
      X_dgps[[dgp]]$forecast_results$temp_name[[sim_run]] <- forecast_results
      
      # go to next sim_run:
      data_start <- data_start + ts_length
      data_end <- data_end + ts_length
      
    }
    
    # rename forecast results list
    new_name <- "hist_avg"
    
    forecast_num <- length(X_dgps[[dgp]]$forecast_results)
    names(X_dgps[[dgp]]$forecast_results) <- c(names(X_dgps[[dgp]]$forecast_results)[-forecast_num], 
                                               new_name)
    
  }
  
  return(X_dgps)
    
}


get_MSE <- function(results_array, y_array, method = "EW"){
  
  # creates MSE array out of results array
  
  # save MSE results in array:
  # (r_squared, correl, coef_vector, model, sim_run)
  # dim nums = (5, 3, 2, 13, 100)
  # get dim of results array and lengths
  res_dims <- dim(results_array)
  r_squared_n <- res_dims[1]
  correls_n <- res_dims[2]
  coef_vector_n <- res_dims[3]
  model_n <- res_dims[4]
  sim_run_n <- res_dims[5]
  # preset array
  if (method == "OW"){
    estim_rat_n <- res_dims[6]
  }
  
  MSE_array <- array(data = NA, dim = res_dims[-length(res_dims)])
  
  
  # calc MSE for each model and save them in MSE array
    
  for (r_squared in c(1:r_squared_n)){
    for (correl in c(1:correls_n)){
      for (coef in c(1:coef_vector_n)){
        for (model in c(1:model_n)){
          for (sim_run in c(1:sim_run_n)){
            
            if (method == "EW"){
              
              MSE_array[r_squared, correl, coef, model, sim_run] <- calc_MSE(predictions = results_array[r_squared, correl, coef, model, sim_run,],
                                                                             targets = y_array[r_squared, correl, coef, sim_run,])
            } else if (method == "OW"){
              
              for (estim_rat in c(1:estim_rat_n)){
                
                MSE_array[r_squared, correl, coef, model, sim_run, estim_rat] <- calc_MSE(predictions = results_array[r_squared, correl, coef, model, sim_run, estim_rat,],
                                                                                          targets = y_array[r_squared, correl, coef, sim_run,])
              }
              
            }
          }
        }
      }
    }
  }
    
  
  return(MSE_array)
}

MSE_normalize <- function(MSE_array, OLS_index = 11){
  
  dimnums <- dim(MSE_array)
  
  # calc MSE for each model and save them in MSE array
  for (r_squared in c(1:dimnums[1])){
    for (correl in c(1:dimnums[2])){
      for (coef in c(1:dimnums[3])){
        for (model in c(1:dimnums[4])){
          for (sim_run in c(1:dimnums[5])){
            MSE_array[r_squared, correl, coef, model, sim_run] <- MSE_array[r_squared, correl, coef, model, sim_run] / MSE_array[r_squared, correl, coef, OLS_index, sim_run]
          }
        }
      }
    }
  }
  
  return(MSE_array)
}






#### functions for combination weight estimation


rolling_ow_estim <- function(x_data, coefs, targets, forecast_start = 71, total_length = 300, weight_estim_ratio = 0.3, non_zero_weights = 0){
  
  # coef start should match x_data start
  # coefs should be 279 x 10 matrix
 
  # initialize values to calculate indexes later
  init_train <- total_length - nrow(coefs)
  
  # preset matrix to save optimal weight estimates
  ow_vector_length <- total_length - forecast_start + 1
  ow_matrix <- matrix(data = NA, nrow = length(coefs[1,]), ncol = ow_vector_length)
  
  for (i in c(1:ow_vector_length)){
    
    x_estim_wind_end <- forecast_start -2 + i
    x_estim_wind_start <- ceiling(x_estim_wind_end*(1-weight_estim_ratio)) + 1
    
    # get oos forecasts
    forecasts <- oos_forecast(t(x_data[x_estim_wind_start:x_estim_wind_end,]), as.matrix(coefs[x_estim_wind_start - init_train -1,]))
    
    # estimate optimal weights
    ow_matrix[,i] <- ow_estim(forecasts, targets[x_estim_wind_start:x_estim_wind_end], non_zero_weights = non_zero_weights)
    
  }
  
  return(ow_matrix)
  
}

oos_forecast <- function(x_data, coefs){
  
  # forecast with x_data as predictor values and coefs as predictor betas
  
  # inputs:
  # x_data is a mastrix of size pxn, n is forecast length and p is num of predictors
  # coefs is a vector of length p
  estims <- matrix(data = NA, nrow = nrow(x_data), ncol = ncol(x_data))
  
  for (data_point in c(1:ncol(x_data))){
    estims[,data_point] <- coefs*x_data[,data_point]
  }
  
  return(estims)
  
}


ow_estim <- function(forecasts, targets, non_zero_weights = 0){
  
  # estimates optimal weights of the forecasts on targets.
  #  weights sum to one
  
  # transpose forecasts
  forecasts <- t(forecasts)
  
  # calc num of vars
  num_vars <- ncol(forecasts)
  
  
  # estimate weights. Last column of forecasts is substracted so that weights sum to one
  lin_model <- lm(formula = (targets - forecasts[,num_vars]) ~.+0, data = as.data.frame(forecasts[, - num_vars] - forecasts[,num_vars]))
  
  # create weights vector from results
  non_last_weights <- as.vector(lin_model$coefficients)
  last_weight <- 1 - sum(non_last_weights)
  weights <- c(non_last_weights, last_weight)
  
  # set negative weights at zero if non_zero_weights = 1
  if (non_zero_weights == 1){
    
    # loop over weights to check if they are negative
    for (i in c(1:length(weights))){
      
      if (weights[i] < 0){
        # set weight to zero
        weights[i] <- 0
      }
    }
    
    # resize weights so that they equal 1
    sum_of_weights <- sum(weights)
    if (sum_of_weights != 1){
      weights <- weights/sum_of_weights
    }
  }
  
  # return weights
  return(weights)
  
}


ow_estim_dgps <- function(X_dgps, forecast_start, weight_estim_ratio_vector,  non_zero_weights = 0){
  
  tic()
  ts_length <- X_dgps[[1]][[1]][[5]]
  num_sim_runs <- X_dgps[[1]][[1]][[6]]
  num_dgps <- length(X_dgps)
  num_models <- length(X_dgps[[1]][[3]])
  n_vars <- length(X_dgps[[1]][[1]][[1]])
  
  num_model_atr <- length(X_dgps[[1]][[3]][[1]][[1]])
  
  
  
  for (dgp_num in c(1:num_dgps)){
    
    
    for (model_num in c(1:num_models)){
      
      start_index <- 1
      
      for (sim_run in c(1:num_sim_runs)){
        
        ow_list <- list()
        i <- 1
        
        for (estim_rat in weight_estim_ratio_vector){
          
          ow_preds <- rolling_ow_estim(x_data = X_dgps[[dgp_num]][[2]][[1]][start_index:(start_index+ts_length -1),1:n_vars], 
                                       coefs = apply(X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]][[4]][,2:(n_vars+1),], c(1,2), sum),
                                       targets = X_dgps[[dgp_num]][[2]][[2]][1:ts_length], total_length = ts_length, 
                                       forecast_start = forecast_start, weight_estim_ratio = estim_rat,
                                       non_zero_weights = non_zero_weights)
          
          ow_list[[i]] <- ow_preds
          i <- i + 1
          
        }
        
        # name optimal weights list with weight est ratios
        names(ow_list) <- as.character(weight_estim_ratio_vector)
        
        # save optimal weights
        X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]][[num_model_atr + 1]] <- ow_list
        
        # rename saved opt weight list based on weight estimation method
        if (non_zero_weights == 1){
          names(X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]]) <-  c(names(X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]])[1:num_model_atr], "non_negative_OWs")
        } else if (non_zero_weights == 0){
          names(X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]]) <-  c(names(X_dgps[[dgp_num]][[3]][[model_num]][[sim_run]])[1:num_model_atr], "unrestricted_OWs")
        }
        
        start_index <- start_index + ts_length
        
          
      }
    }
  }
  toc()
  return(X_dgps)
  
}




get_results <- function(X_dgps, r_squareds = c(0.01,0.025,0.05,0.1,0.25),
                        correls = c(0.1,0.5,0.8), 
                        coefs = list("full" = rep(1,times =10), 
                                     "sparse" = c(rep(1,times =5), rep(0, times = 5))), 
                        method = "fc", fc_length = 230, ow_type_index = 5){
  # save results in array:
  # (r_squared, correl, coef_vector, model, sim_run, estim_ratio)
  # dim nums = (5, 3, 1, 13, 100, 3)
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  # fc_length <- length(X_dgps[[1]][[3]][[1]][[1]][[1]])
  estimated_fcs_length <- length(X_dgps[[1]][[3]][[1]][[1]][[1]])
  fc_start <- estimated_fcs_length - fc_length + 1
  if (method == "OW"){
    estim_ratio_dim <- length(X_dgps[[1]][[3]][[1]][[1]][[5]])
  }
  
  
  if (method == "fc"){
    
    model_num <- length(X_dgps[[1]][[3]])
    dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, sim_num, fc_length)
    
  } else if (method == "targets"){
    
    dimnums <- c(r_squared_dim, correl_dim, coefs_dim, sim_num, fc_length)
    
    len <- X_dgps[[1]][[1]][[5]]
  } else if (method == "OW"){
    
    model_num <- length(X_dgps[[1]][[3]])
    dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, sim_num, estim_ratio_dim, fc_length)
  }
  
  result_array <- array(data = NA, dim = dimnums)
  
  
  
  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    
    if (method == "fc"){
      # create array of forecasts
      
      # iterate over models
      for (model_index in c(1:model_num)){
        # iterate over sim runs
        for (sim_run_index in c(1:sim_runs)){
          # iterate over data point forecasts
          
          for (fc_index in  c(fc_start:estimated_fcs_length)){
            
            # save data point forecast
            result_array[r_squred_index, correl_index, coef_index, model_index, sim_run_index, fc_index - fc_start + 1] <- 
              dgp[[3]][[model_index]][[sim_run_index]][[1]][[fc_index]]
            
          }
          
        }
      }
    } else if (method == "targets"){
      # create array of y data
      data_start <- len - fc_length
      
      # iterate over sim runs
      for (sim_run_index in c(1:sim_runs)){
        # iterate over data points
        for (fc_index in c(1:fc_length)){
          # save y value
          result_array[r_squred_index, correl_index, coef_index, sim_run_index, fc_index] <- 
            dgp[[2]][[2]][[fc_index+data_start]]
        }
        
        # jump to the point in the data where next forecast starts
        data_start <- data_start + len
      }
      
      
    } else if (method == "OW"){
      # create array of forecasts
      
      # iterate over models
      for (model_index in c(1:model_num)){
        # iterate over sim runs
        for (sim_run_index in c(1:sim_runs)){
          
          # iterate over estim ratios
          for (estim_rat_index in c(1:estim_ratio_dim)){
            
            # iterate over data point forecasts
            for (fc_index in  c(fc_start:estimated_fcs_length)){
              
              # get OWs
              OWs <- dgp[[3]][[model_index]][[sim_run_index]][[ow_type_index]][[estim_rat_index]][,fc_index - fc_start + 1]
              
              # get indiv forecasts
              ind_fcs <- dgp[[3]][[model_index]][[sim_run_index]][[3]][,fc_index]
              
              # calculate forecast with optimal weights
              OW_fc <- sum(OWs * ind_fcs)
              
              # save data point forecast
              result_array[r_squred_index, correl_index, coef_index, model_index, sim_run_index, 
                           estim_rat_index, fc_index - fc_start + 1] <- OW_fc
              
            }
            
          }
          
          
        }
      }
    }
    
  }
  
  return(result_array)
}


get_ow_r_squareds <- function(X_dgps, r_squareds = c(0.01,0.025,0.05,0.1,0.25),
                        correls = c(0.1,0.5,0.8), 
                        coefs = list("full" = rep(1,times =10), 
                                     "sparse" = c(rep(1,times =5), rep(0, times = 5)))){
  # save results in array:
  # (r_squared, correl, coef_vector, model, sim_run, estim_ratio)
  # dim nums = (5, 3, 1, 13, 100, 3)
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  num_of_vars <- X_dgps[[1]][[1]][[4]]
  
    
  model_num <- length(X_dgps[[1]][[3]])
  dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num)
  
  r_squared_array <- array(data = 0, dim = dimnums)
  
  
  
  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    

    # iterate over models
    for (model_index in c(1:model_num)){
      # iterate over sim runs
      for (sim_run_index in c(1:sim_num)){
        # iterate over data point forecasts
          
        x_data <- t(dgp[[3]][[model_index]][[sim_run_index]][[3]])
          
        x_data <- x_data[,-num_of_vars] - x_data[,num_of_vars]
          
        x_data <- as.data.frame(x_data)
        
        lin_model <- lm(formula = V1 ~.+0, data = x_data)  
          
        r_squared_array[r_squred_index, correl_index, coef_index, model_index] <- 
          r_squared_array[r_squred_index, correl_index, coef_index, model_index] + summary(lin_model)$r.squared
            
      }
    }
  }
  
  
  r_squared_array <- r_squared_array / sim_num
  
  return(r_squared_array)
}



create_hist_array <- function(X_dgps, coefs = c(1,1,1), ow_index = 5, sim_num_max = 100, num_var_max = 10){
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  if (sim_num > sim_num_max){
    sim_num <- sim_num_max
  }
  num_of_vars <- X_dgps[[1]][[1]][[4]]
  if (num_of_vars > num_var_max){
    num_of_vars <- num_var_max
  }
  ow_len <- length(X_dgps[[1]][[3]][[1]][[1]][[ow_index]][[1]][1,])
  hist_len <- ow_len*sim_num
  weight_estim_rat_dim <- 3
  
  
  model_num <- length(X_dgps[[1]][[3]])
  dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, num_of_vars, hist_len, weight_estim_rat_dim)
  
  hist_array <- array(data = NA, dim = dimnums)
  

  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    
    
    # iterate over models
    for (model_index in c(1:model_num)){
      # iterate over sim runs
      for (sim_run_index in c(1:sim_num)){
        # iterate over data point forecasts
        
        for (var in c(1:num_of_vars)){
          
          for (weight_estim_rat_ind in c(1:weight_estim_rat_dim)){
            hist_array[r_squred_index, correl_index, coef_index, model_index, var, c(((sim_run_index-1)*ow_len +1):(sim_run_index*ow_len)),weight_estim_rat_ind] <- 
              dgp[[3]][[model_index]][[sim_run_index]][[ow_index]][[weight_estim_rat_ind]][var,]
            
          }
          
        }
      }
    }
  }
  
  return(hist_array)
  
  
}



correl_calc <- function(X_dgps, coefs = c(1,1,1)){
  
  # save results in array:
  # (r_squared, correl, coef_vector, model, sim_run, estim_ratio)
  # dim nums = (5, 3, 1, 13, 100, 3)
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  num_of_vars <- X_dgps[[1]][[1]][[4]]
  
  
  model_num <- length(X_dgps[[1]][[3]]) -1
  dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, 7)
  
  correl_array <- array(data = 0, dim = dimnums)
  #browser()
  
  
  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    
    
    # iterate over models
    for (model_index in c(1:model_num)){
      # iterate over sim runs
      for (sim_run_index in c(1:sim_num)){
        # iterate over data point forecasts
        
        x_data <- t(dgp[[3]][[model_index]][[sim_run_index]][[3]])
        
        correl <- cor(x_data)
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 1] <- 
          correl[1,2] + correl_array[r_squred_index,correl_index,coef_index, model_index, 1]
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 2] <- 
          correl[1,10] + correl_array[r_squred_index,correl_index,coef_index, model_index, 2]
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 3] <- 
          correl[2,10] + correl_array[r_squred_index,correl_index,coef_index, model_index, 3]
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 4] <- 
          sum(correl[upper.tri(correl)])/length(correl[upper.tri(correl)]) + correl_array[r_squred_index,correl_index,coef_index, model_index, 4]
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 5] <- 
          sum(abs(correl[upper.tri(correl)]))/length(correl[upper.tri(correl)]) + correl_array[r_squred_index,correl_index,coef_index, model_index, 5]
        
        
        x_data <- x_data[,-num_of_vars] - x_data[,num_of_vars]
        
        correl <- cor(x_data)
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 6] <- 
          sum(correl[upper.tri(correl)])/length(correl[upper.tri(correl)]) + correl_array[r_squred_index,correl_index,coef_index, model_index, 6]
        
        correl_array[r_squred_index,correl_index,coef_index, model_index, 7] <-
          sum(abs(correl[upper.tri(correl)]))/length(correl[upper.tri(correl)]) + correl_array[r_squred_index,correl_index,coef_index, model_index, 7]
        
      }
    }
  }
  
  
  correl_array <- correl_array / sim_num
  
  return(correl_array)
}





get_cor_mat <- function(X_dgps, coefs = c(1,1,1), method = "normal"){
  
  # save results in array:
  # (r_squared, correl, coef_vector, model, sim_run, estim_ratio)
  # dim nums = (5, 3, 1, 13, 100, 3)
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  num_of_vars <- X_dgps[[1]][[1]][[4]]
  
  
  model_num <- length(X_dgps[[1]][[3]])
  if (method == "normal"){
    cor_n <- num_of_vars^2
  } else if (method == "minus"){
    cor_n <- (num_of_vars-1)^2
  }
  
  dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, cor_n)
  
  cor_mat_array <- array(data = 0, dim = dimnums)
  
  
  
  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    
    
    # iterate over models
    for (model_index in c(1:model_num)){
      # iterate over sim runs
      for (sim_run_index in c(1:sim_num)){
        # iterate over data point forecasts
        
        x_data <- t(dgp[[3]][[model_index]][[sim_run_index]][[3]])
        
        if (method == "normal"){
          correl <- cor(x_data)
          
          for (i in c(1:cor_n)){
            cor_mat_array[r_squred_index,correl_index,coef_index, model_index, i] <- 
              correl[i] + cor_mat_array[r_squred_index,correl_index,coef_index, model_index, i]
          }
        } else if (method == "minus"){
            
          correl <- cor(x_data[,-num_of_vars] - x_data[, num_of_vars])
          
          for (i in c(1:cor_n)){
            cor_mat_array[r_squred_index,correl_index,coef_index, model_index, i] <- 
              correl[i] + cor_mat_array[r_squred_index,correl_index,coef_index, model_index, i]
          }
        }
          
        }
        
        
      }
    }
  
  
  cor_mat_array <- cor_mat_array / sim_num
  
  return(cor_mat_array)
}




validate_lambda <- function(X_dgps, valid_start, valid_length, lambda_num = 10, method = NA){
  
  # get length vals from X_dgps
  dgp_num <- length(X_dgps)
  ts_length <- X_dgps[[1]][[1]][[5]]
  sim_runs_num <- X_dgps[[1]][[1]][[6]]
  
  forecast_length <- length(X_dgps[[1]][[3]][[1]][[1]][[1]])
  valid_start_corrigated <- valid_start - (ts_length - forecast_length)
  
  # preset results matrix
  validated_matrix <- matrix(data = NA, nrow = 3, ncol = ts_length - (valid_start + valid_length -1))
  
  for (dgp_n in c(1:dgp_num)){
    
    # for indexing y data later
    y_drift <- 0
    
    for (sim_run in c(1:sim_runs_num)){
      
      # for indexin validated_matrix later
      data_point_drift <- 0
      
      for (data_point in c((valid_start+valid_length):ts_length)){
        
        # preset MSE vector
        MSE_val <- vector(length = lambda_num)
        # preset start and end indices
        fc_start <- valid_start_corrigated+ data_point_drift
        fc_end <- valid_start_corrigated+valid_length-1 + data_point_drift
        y_start <- valid_start + data_point_drift + y_drift*ts_length
        y_end <- valid_start + data_point_drift + valid_length-1 + y_drift*ts_length
        
        if (method == "expanding"){
          fc_start <- fc_start - data_point_drift
          y_start <- y_start - data_point_drift
        }
      
        for (lambda_n in c(1:lambda_num)){
        
          # get forecast and target data
          forecast_data <- X_dgps[[dgp_n]][[3]][[lambda_n]][[sim_run]][[1]][fc_start:fc_end]
          y_data <- X_dgps[[dgp_n]][[2]][[2]][y_start:y_end]
          
          # calc MSE on val data
          MSE_val[lambda_n] <- calc_MSE(forecast_data, y_data)
        }
        
        # save: data point index, optimal lambda, optimal forecast
        validated_matrix[1,data_point_drift+1] <- data_point
        validated_matrix[2,data_point_drift+1] <- which.min(MSE_val)/lambda_num
        validated_matrix[3,data_point_drift+1] <- X_dgps[[dgp_n]][[3]][[which.min(MSE_val)]][[sim_run]][[1]][[data_point - (ts_length - forecast_length)]]
        
        # drift data point index
        data_point_drift <- data_point_drift + 1
      }
      
      # drift y indexing
      y_drift <- y_drift + 1
      
      X_dgps[[dgp_n]]$validated_lambda[[sim_run]] <- validated_matrix 
    }
    
    # rename new part as "validated lambda"
    #names(X_dgps[[dgp_n]]) <- c(names(X_dgps[[dgp_n]])[-4], "validated_lambda")
  }
  
    
  return(X_dgps)
  
}

get_validated_res <- function(X_dgps, correls = c(0,0.4,0.8), 
                              r_squareds = c(0.01,0.025,0.05,0.1,0.25),
                              coefs = rep(1, times = 10)){
  
  
  
  r_squared_dim <- length(r_squareds)
  correl_dim <-  length(correls)
  coefs_dim <- 1
  sim_num <- X_dgps[[1]][[1]][[6]]
  fc_length <- ncol(X_dgps[[1]][[4]][[1]])
  #estimated_fcs_length <- length(X_dgps[[1]][[3]][[1]][[1]][[1]])
  #fc_start <- estimated_fcs_length - fc_length + 1
  model_num <- 1

  dimnums <- c(r_squared_dim, correl_dim, coefs_dim, model_num, sim_num, fc_length)


  result_array <- array(data = NA, dim = dimnums)

  for (dgp in X_dgps){
    
    # get dgp params
    r_squared <- dgp$dgp_param$r_squared
    correl <- dgp$dgp_param$correl
    coef <- dgp$dgp_param$coefs
    
    # find dgp param indices
    r_squred_index <- match(r_squared,r_squareds)
    correl_index <- match(correl, correls)
    coef_index <- 1
    
      # create array of forecasts
      
        # iterate over sim runs
        for (sim_run_index in c(1:sim_runs)){
          # iterate over data point forecasts
          
          for (fc_index in  c(1:fc_length)){
            
            # save data point forecast
            result_array[r_squred_index, correl_index, coef_index, model_num, sim_run_index, fc_index] <- 
              dgp[[4]][[sim_run_index]][3,fc_index]
            
          }
          
        }
      
    } 
    
    return(result_array)
    
}


calc_bias_var_covar <- function(predictions_list, targets){
  #browser()
  # preset vectors for storage
  len <- length(predictions_list)
  bias_squared <- vector(length = len)
  variance <- vector(length = len)
  covariance <- vector(length = len*(len-1)/2)
  i <- 1
  k <- 1
  
  for (predictor_model in predictions_list){
    
    bias_squared[i] <- (mean(targets - predictor_model))^2
    variance[i] <- sum((rep(mean(predictor_model), times = length(predictor_model)) - predictor_model)^2)/(length(predictor_model)-1)
    
    for (j in c(1:len)){
      if (j > i){
        covariance[k] <- sum((rep(mean(predictor_model), times =length(predictor_model)) - predictor_model) * (rep(mean(predictions_list[[j]]), times = length(predictor_model)) - predictions_list[[j]]))/(length(predictor_model)-1)

          #cov(predictor_model, predictions_list[[j]])
        k <- k + 1
      }
    }
  }
  
  return_vector <- vector(length = 3)
  
  return_vector[1] <- mean(bias_squared)
  return_vector[2] <- mean(variance) / len
  return_vector[3] <- mean(covariance) * (1 - 1/len)
  
  return(return_vector)
}









  
  
asset_allocation <- function(actual_excess_rets, risk_frees, ret_forecast, volatility_forecast,
                             gamma, c_bp){
  
  
  
  # Input
  #
  # actual        = T-vector of actual excess returns
  # risk_free     = T-vector of risk-free rates
  # forecast      = T-vector of excess return forecasts
  # volatility_FC = T-vector of volatility forecasts
  # gamma_MV      = risk aversion parameter
  # c_bp          = transaction cost in basis points
  #
  # Output
  #
  # avg_utility       = average utility
  # SR                = Sharpe ratio
  # weight_risky      = T-vector of risky portfolio weights
  # cumulative_return = T-vector of cumulative portfolio returns
  # avg_turnover      = average turnover
  #browser()
  c <- c_bp/10000
  inv_length <- length(actual_excess_rets)
  weight_risky <- vector(length = inv_length)
  turnover <- vector(length = inv_length)
  return_portfolio <- vector(length = inv_length)
  
  for (time in c(1:inv_length)){
    weight_risky[time] <- (1/gamma)*(ret_forecast[time]/volatility_forecast[time])
    
    if (weight_risky[time] < 0){
      weight_risky[time] <- 0
    } else if (weight_risky[time] > 1.5){
      weight_risky[time] <- 1.5
    }
  }
  
  for (time in c(1:inv_length)){
    
    if (time < inv_length){
      wealth_total_end_t <- 1 + risk_frees[time] + weight_risky[time]*actual_excess_rets[time]
      wealth_risky_end_t <- weight_risky[time]*(1+risk_frees[time]+actual_excess_rets[time])
      target_risky_end_t <- weight_risky[time+1]*wealth_total_end_t
      turnover[time] <- abs(target_risky_end_t - wealth_risky_end_t)/wealth_total_end_t
      TC_t = c*turnover[time]
      return_portfolio[time] <- (1+risk_frees[time] + weight_risky[time]*actual_excess_rets[time])/(1-c*turnover[time]) - 1
    } else {
      return_portfolio[time] <- risk_frees[time] + weight_risky[time]*actual_excess_rets[time]
    }

  }
  
  return_list <- list()
  return_list$avg_utility <- mean(return_portfolio) - 0.5*gamma*(sd(return_portfolio))^2
  return_list$excess_return <- return_portfolio - risk_frees
  return_list$SR <- mean(return_list$excess_return)/sd(return_list$excess_return)
  return_list$cumulative_return <- cumprod(1+return_portfolio)
  return_list$avg_turnover <- mean(turnover)
  
  return(return_list)

}

volatility_forecast <- function(stock_return_data, estimation_window = 40){
  
  # forecasts volatility of stock return form stock return data
  # with estimation window length of 'estimation_window'
  
  vol_forecast <- vector(length = length(stock_return_data)-estimation_window)
  
  for (estim_index in c(1:(length(stock_return_data)-estimation_window))){
    
    vol_forecast[estim_index] <- var(stock_return_data[estim_index:(estim_index+estimation_window)])
  }
  
  return(vol_forecast)
}












