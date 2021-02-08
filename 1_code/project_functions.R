
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
  x <- mvtnorm::rmvnorm(ts_length, mean = rep(0, times = num_of_vars), 
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
  
  y <- rowSums(t(resized_coefs * t(x))) + rnorm(ts_length, mean = 0, sd = 1)
  
  # store generated time series in storage matrix X
  X[[dgp_index]]$sim_data$x <- x
  X[[dgp_index]]$sim_data$y <- y
  
  # store new params in dgp_param
  X[[dgp_index]]$dgp_param$resized_coefs <- resized_coefs
  X[[dgp_index]]$dgp_param$covar_matrix <- covar_matrix
  
  return(X)
  
}

### próba
X <- list()
X <- def_dgp_params(X, coefs = c(1,1,0,0), r_squared = 0.9, correl = 0.3, num_of_vars = 4, len = 500)
X <- simulate_dgp(X,1)

