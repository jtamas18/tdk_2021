
# set working directory if necessary first
setwd("C:/Users/user/Desktop/git/tdk_2021")

source("1_code/project_functions.R")


library("readxl")
library("ggplot2")
library("gridExtra")

#### DATA PREPARATION ####


original_data <- read_excel(path = "0_data/data_goyal.xlsx", sheet = 2, col_types = "numeric")

# drop early data
original_data <- original_data[-c(1:302),]

# lag index
original_data$lagged_index <- NA
for(index in c(2:length(original_data$Index))){
  original_data$lagged_index[index] <- original_data$Index[index-1]
}

# lag risk-free rate 
original_data$lagged_rf <- NA
for(index in c(2:length(original_data$Rfree))){
  original_data$lagged_rf[index] <- original_data$Rfree[index-1]
}


# calculate equity risk premium
original_data$equity_premium <- original_data$CRSP_SPvw - original_data$lagged_rf
original_data$log_equity_premium <- log(1 + original_data$CRSP_SPvw) - log(1 + original_data$lagged_rf)

#original_data$IndexDiv <- original_data$Index + original_data$D12
#original_data$return_stocks <- log(original_data$IndexDiv) - log(original_data$lagged_index)
#original_data$equity_premium <- original_data$return_stocks - log(original_data$Rfree + 1)

# lag equity premium
original_data$next_log_equity_prem <- NA
for(index in c(2:length(original_data$next_log_equity_prem))){
  original_data$next_log_equity_prem[index] <- original_data$log_equity_premium[index+1]
}

original_data$next_equity_prem <- NA
for(index in c(2:length(original_data$next_equity_prem))){
  original_data$next_equity_prem[index] <- original_data$equity_premium[index+1]
}



# preset data table and add index and dividends

data <- data.frame(equity_premium = original_data$next_log_equity_prem)
#data$non_log_equity_premium <- original_data$next_equity_prem

# add date as rownames
rownames(data) <- original_data$yyyyq

# dividend price ratio
data$dp <- log(original_data$D12) - log(original_data$Index) 



# divindend yield
data$dy <- log(original_data$D12) - log(original_data$lagged_index)

# earning price ratio
data$ep <- log(original_data$E12) - log(original_data$Index)

# dividend payout (dividend to earnings) ratio
data$de <- log(original_data$D12) - log(original_data$E12)

# stock variance (svar)
data$svar <- original_data$svar

# book to market ratio
data$bm <- original_data$`b/m`

# net equity expnasion
data$ntis <- original_data$ntis

# treasury bills
data$tbl <- original_data$tbl

# long term yield
data$lty <- original_data$lty

# long term rate of return
data$ltr <- original_data$ltr 

# term spread
data$tms <- original_data$lty - original_data$tbl

# default yield spread
data$dfy <- original_data$BAA - original_data$AAA

# default return spread
data$dfr <- original_data$corpr - original_data$ltr

# inflation
  # lag inflation
  original_data$lagged_infl <- NA
  for(index in c(2:length(original_data$infl))){
    original_data$lagged_infl[index] <- original_data$infl[index-1]
  }

data$infl <- original_data$lagged_infl

# investment to capital ratio
data$ik <- original_data$ik

# save risk frees
risk_frees <- original_data$Rfree
risk_frees <- risk_frees[-c(1,2,length(risk_frees))]

# save initial stock returns
actual_stock_returns <- original_data$CRSP_SPvw
actual_stock_returns <- actual_stock_returns[-c(1,2,length(actual_stock_returns))]

# drop observations before 1947Q1
data <- data[-c(1,2),]
data <- data[-nrow(data),] # equity premium is NA here

rm(original_data)

#### MODEL FITTING ####

# get data
eq_prem <- data$equity_premium
predictor_data <- as.matrix(data[,-1])

selector_mat_list <- choose_models(predictor_data, method = "univariate")

init_train <- 72

# fit models

univariate_models <- list()

univariate_models$trad <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                               selector_mat_list = selector_mat_list, init_train = init_train,
                               incl_cons = 1, method = "traditional", l_rate = 0.18, n_iter = 10000)
univariate_models$ncl_0.1 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.1, l_rate = 0.10, n_iter = 150000,
                                         OLS_coefs = univariate_models$trad$ind_coef)


univariate_models$ncl_0.2 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.2, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.1$ind_coef)


univariate_models$ncl_0.3 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.3, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.2$ind_coef)

univariate_models$ncl_0.4 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.4, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.3$ind_coef)

univariate_models$ncl_0.5 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.5, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.4$ind_coef)


univariate_models$ncl_0.6 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.6, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.5$ind_coef)

univariate_models$ncl_0.7 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.7, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.6$ind_coef)


univariate_models$ncl_0.8 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.8, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.7$ind_coef)

univariate_models$ncl_0.9 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 0.9, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.8$ind_coef)

univariate_models$ncl_1.0 <- forecast_fc(y_data = eq_prem, x_data = predictor_data, 
                                         selector_mat_list = selector_mat_list, init_train = init_train,
                                         incl_cons = 1, method = "ncl", lambda = 1.0, l_rate = 0.10, n_iter = 100000,
                                         OLS_coefs = univariate_models$ncl_0.9$ind_coef)



hist_avg_forecasts <- forecast_hist_avg(y_data = eq_prem, init_train = init_train)







## EVALUATE MODELS ON OOS DATA

eq_prem_oos <- eq_prem[-c(1:init_train)]

# calc MSEs

num_of_models <- length(univariate_models)

MSEs <- vector(length = (num_of_models+1))

for (i in c(1:num_of_models)){
  MSEs[i] <- calc_MSE(predictions = univariate_models[[i]][[1]], targets = eq_prem_oos)
}

MSEs[num_of_models+1] <- calc_MSE(predictions = hist_avg_forecasts, targets = eq_prem_oos)

### calculate oos R^2-s

R_squared <- vector(length = num_of_models)

for (i in c(1:num_of_models)){
  R_squared[i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts,
                                 predictions = as.vector(univariate_models[[i]][[1]]))
}


## clark-west test

cw_p_vals <- vector(length = num_of_models)

for (i in c(1:num_of_models)){
  cw_p_vals[i] <- clark_west_test(targets = eq_prem_oos, predictions = univariate_models[[i]][[1]],
                               hist_avg_pred = hist_avg_forecasts)[2]
  
}




## calc Campbell restrictions

# 

positive_coefs <- vector(length = 15)
for (i in c(1:15)){
  if (sum(sign(univariate_models$trad$fc_coef[,i+1])) > 0){
    positive_coefs[i] <- 1
  } else {
    positive_coefs[i] <- -1
  }
}

univar_campbell <- univariate_models


for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    for (coef in c(1:15)){
      
      # test if coef is not of the expected sign
      if (sign(univariate_models[[model_num]][[2]][fc,coef+1]) != positive_coefs[coef]){
        # set to zero if not epexted sign
        univar_campbell[[model_num]][[2]][fc,coef+1] <- 0
      }
    }
  }
}

# recalc forecasts
oos_pred_data <- predictor_data[-c(1:init_train),]

for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    # calc forecasts with new coefs
    univar_campbell[[model_num]][[1]][fc] <- 
      sum(oos_pred_data[fc,]*univar_campbell[[model_num]][[2]][fc,2:16]) +
      univar_campbell[[model_num]][[2]][fc,1]
  }
}


# set negative forecasts to zero

for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    # check if forecast is negative
    # set to zero if it is
    if (univar_campbell[[model_num]][[1]][fc] < 0){
      univar_campbell[[model_num]][[1]][fc] <- 0
    }
  }
}

## evaluate Campbell restricted models

eq_prem_oos <- eq_prem[-c(1:init_train)]

# calc MSEs

num_of_models <- length(univar_campbell)

MSEs_Campbell <- vector(length = (num_of_models+1))

for (i in c(1:num_of_models)){
  MSEs_Campbell[i] <- calc_MSE(predictions = univar_campbell[[i]][[1]], targets = eq_prem_oos)
}

MSEs_Campbell[num_of_models+1] <- calc_MSE(predictions = hist_avg_forecasts, targets = eq_prem_oos)



## plot nonrestricted and Campbell restricted

par(mfrow=c(1,3))

plot(MSEs)
plot(MSEs_Campbell)
plot(c(MSEs,MSEs_Campbell))
plot(MSEs_nonneg)
plot(c(MSEs,MSEs_nonneg,MSEs_coef[c(1:10)]))

# nonneq equity premium restriction


univar_nonneg <- univariate_models


# set negative forecasts to zero

for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    # check if forecast is negative
    # set to zero if it is
    if (univar_nonneg[[model_num]][[1]][fc] < 0){
      univar_nonneg[[model_num]][[1]][fc] <- 0
    }
  }
}



## evaluate nonneg restricted models

eq_prem_oos <- eq_prem[-c(1:init_train)]

# calc MSEs

num_of_models <- length(univar_nonneg)

MSEs_nonneg <- vector(length = (num_of_models+1))

for (i in c(1:num_of_models)){
  MSEs_nonneg[i] <- calc_MSE(predictions = univar_nonneg[[i]][[1]], targets = eq_prem_oos)
}

#MSEs_nonneg[num_of_models+1] <- calc_MSE(predictions = hist_avg_forecasts, targets = eq_prem_oos)


plot(c(MSEs,MSEs_Campbell,MSEs_nonneg))





# coef restriction


positive_coefs <- vector(length = 15)
for (i in c(1:15)){
  if (sum(sign(univariate_models$trad$fc_coef[,i+1])) > 0){
    positive_coefs[i] <- 1
  } else {
    positive_coefs[i] <- -1
  }
}

univar_coef <- univariate_models


for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    for (coef in c(1:15)){
      
      # test if coef is not of the expected sign
      if (sign(univariate_models[[model_num]][[2]][fc,coef+1]) != positive_coefs[coef]){
        # set to zero if not epexted sign
        univar_coef[[model_num]][[2]][fc,coef+1] <- 0
      }
    }
  }
}

# recalc forecasts
oos_pred_data <- predictor_data[-c(1:init_train),]

for (model_num in c(1:length(univariate_models))){
  
  for (fc in c(1:length(univariate_models$trad$fc_forecast))){
    
    # calc forecasts with new coefs
    univar_coef[[model_num]][[1]][fc] <- 
      sum(oos_pred_data[fc,]*univar_coef[[model_num]][[2]][fc,2:16]) +
      univar_coef[[model_num]][[2]][fc,1]
  }
}




## evaluate coef restricted models

eq_prem_oos <- eq_prem[-c(1:init_train)]

# calc MSEs

num_of_models <- length(univar_coef)

MSEs_coef <- vector(length = (num_of_models+1))

for (i in c(1:num_of_models)){
  MSEs_coef[i] <- calc_MSE(predictions = univar_coef[[i]][[1]], targets = eq_prem_oos)
}

MSEs_coef[num_of_models+1] <- calc_MSE(predictions = hist_avg_forecasts, targets = eq_prem_oos)

plot(c(MSEs,MSEs_Campbell,MSEs_nonneg,MSEs_coef[c(1:10,12)]))






# Rapach (2010) test

## evaluate models on oos sample

eq_prem_oos <- eq_prem[c((init_train+1):236)]

# calc MSEs

num_of_models <- length(univariate_models)

MSEs <- vector(length = (num_of_models+1))

for (i in c(1:num_of_models)){
  MSEs[i] <- calc_MSE(predictions = univariate_models[[i]][[1]][c(1:(236-init_train))], targets = eq_prem_oos)
}

MSEs[num_of_models+1] <- calc_MSE(predictions = hist_avg_forecasts[c(1:(236-init_train))], targets = eq_prem_oos)

### calculate oos R^2-s

R_squared <- vector(length = num_of_models)

for (i in c(1:num_of_models)){
  R_squared[i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts[c(1:(236-init_train))],
                                 predictions = as.vector(univariate_models[[i]][[1]][c(1:(236-init_train))]))
}

calc_r_squared(targets = eq_prem_oos, predictions = univar_nonneg$ncl_1.0$fc_forecast, hist_avg_pred = hist_avg_forecasts)




















#### VALIDATED LAMBDA ####

init_valid <- 20
fc_length <- length(univariate_models$trad$fc_forecast)
validation_length <- fc_length - init_valid
validated_matrix  <- matrix(data = NA, nrow = 2, 
                           ncol = validation_length)
validated_campbell_matrix <- matrix(data = NA, nrow = 2, 
                                  ncol = validation_length)

validated_nonneg_matrix <- matrix(data = NA, nrow = 2, 
                                  ncol = validation_length)
validated_coef_matrix <- matrix(data = NA, nrow = 2, 
                                  ncol = validation_length)


num_of_models <- length(univariate_models)
MSE_matrix <- matrix(data = NA, nrow = 4, ncol = num_of_models)

eq_prem_oos <- eq_prem[-c(1:init_train)]

i <- 1
for (estim_window_end in c(init_valid:(fc_length-1))){
  
  for (model_index in c(1:num_of_models)){
    
    MSE_matrix[1,model_index] <- calc_MSE(targets = eq_prem_oos[1:estim_window_end], 
                                          predictions = univariate_models[[model_index]][[1]][1:estim_window_end])
    MSE_matrix[2,model_index] <- calc_MSE(targets = eq_prem_oos[1:estim_window_end], 
                                          predictions = univar_campbell[[model_index]][[1]][1:estim_window_end])
    MSE_matrix[3,model_index] <- calc_MSE(targets = eq_prem_oos[1:estim_window_end], 
                                          predictions = univar_nonneg[[model_index]][[1]][1:estim_window_end])
    MSE_matrix[4,model_index] <- calc_MSE(targets = eq_prem_oos[1:estim_window_end], 
                                          predictions = univar_coef[[model_index]][[1]][1:estim_window_end])
    
  }
  
  validated_matrix[1,i] <- (which.min(MSE_matrix[1,c(1:10)]) - 1)/10
  validated_campbell_matrix[1,i] <- (which.min(MSE_matrix[2,]) - 1)/10
  validated_nonneg_matrix[1,i] <- (which.min(MSE_matrix[3,]) - 1)/10
  validated_coef_matrix[1,i] <- (which.min(MSE_matrix[4,]) - 1)/10
  
  validated_matrix[2,i] <- univariate_models[[which.min(MSE_matrix[1,c(1:10)])]][[1]][init_valid+i]
  validated_campbell_matrix[2,i] <- univar_campbell[[which.min(MSE_matrix[2,])]][[1]][init_valid+i]
  validated_nonneg_matrix[2,i] <- univar_nonneg[[which.min(MSE_matrix[3,])]][[1]][init_valid+i]
  validated_coef_matrix[2,i] <- univar_coef[[which.min(MSE_matrix[4,])]][[1]][init_valid+i]
  
  i <- i + 1
  
}


# evaluate validated

validated_r_squared <- vector(length = 4)
validated_r_squared[1] <- calc_r_squared(targets = eq_prem_oos[-c(1:init_valid)], 
                                        predictions = validated_matrix[2,], 
                                        hist_avg_pred = hist_avg_forecasts[-c(1:init_valid)])
validated_r_squared[2] <- calc_r_squared(targets = eq_prem_oos[-c(1:init_valid)], 
                                        predictions = validated_campbell_matrix[2,], 
                                        hist_avg_pred = hist_avg_forecasts[-c(1:init_valid)])
validated_r_squared[3] <- calc_r_squared(targets = eq_prem_oos[-c(1:init_valid)], 
                                        predictions = validated_nonneg_matrix[2,], 
                                        hist_avg_pred = hist_avg_forecasts[-c(1:init_valid)])
validated_r_squared[4] <- calc_r_squared(targets = eq_prem_oos[-c(1:init_valid)], 
                                        predictions = validated_coef_matrix[2,], 
                                        hist_avg_pred = hist_avg_forecasts[-c(1:init_valid)])
validated_r_squared



# rolling R_squareds
rolling_r_squared <- vector(length = validation_length)
for (i in c(1:validation_length)){
  rolling_r_squared[i] <- calc_r_squared(targets = eq_prem_oos[-c(1:(init_valid+i-1))], 
                                         predictions = univar_nonneg$ncl_1.0$fc_forecast[(init_valid+i):length(univariate_models$trad$fc_forecast)], 
                                         hist_avg_pred = hist_avg_forecasts[-c(1:(init_valid+i-1))])*100
}



#### ASSET ALLOCATION

# get stock return predictions

stock_returns <- matrix(data = NA, nrow = 5, ncol = length(univariate_models$trad$fc_forecast))
stock_returns_val <- matrix(data = NA, nrow = 5, ncol = length(validated_matrix[2,]))

stock_returns[1,] <- (exp(univariate_models$trad$fc_forecast)*(1+risk_frees[-c(1:init_train)])-1)
stock_returns[2,] <- exp(univariate_models[[which.min(MSEs[1:11])]][[1]])*(1+risk_frees[-c(1:init_train)])-1
stock_returns[3,] <- exp(univar_campbell[[which.min(MSEs_Campbell[1:11])]][[1]])*(1+risk_frees[-c(1:init_train)])-1
stock_returns[4,] <- exp(univar_nonneg[[which.min(MSEs_nonneg[1:11])]][[1]])*(1+risk_frees[-c(1:init_train)])-1
stock_returns[5,] <- exp(univar_coef[[which.min(MSEs_coef[1:11])]][[1]])*(1+risk_frees[-c(1:init_train)])-1

  
stock_returns_val[1,] <- exp(univariate_models$trad$fc_forecast[-c(1:init_valid)])*(1+risk_frees[-c(1:(init_train+init_valid))])-1
stock_returns_val[2,] <- exp(validated_matrix[2,])*(1+risk_frees[-c(1:(init_train+init_valid))])-1
stock_returns_val[3,] <- exp(validated_campbell_matrix[2,])*(1+risk_frees[-c(1:(init_train+init_valid))])-1
stock_returns_val[4,] <- exp(validated_nonneg_matrix[2,])*(1+risk_frees[-c(1:(init_train+init_valid))])-1
stock_returns_val[5,] <- exp(validated_coef_matrix[2,])*(1+risk_frees[-c(1:(init_train+init_valid))])-1

# predict volatility

vol_forecasts <- volatility_forecast(stock_return_data = actual_stock_returns[c((init_train-40):(length(actual_stock_returns)-1))],
                                     estimation_window = 40)


# calculate asset management

asset_management_sim <- list()

for (i in c(1:length(stock_returns[,1]))){
  asset_management_sim[[i]] <- asset_allocation(actual_excess_rets = actual_stock_returns[-c(1:init_train)] - risk_frees[-c(1:init_train)],
                                                risk_frees = risk_frees[-c(1:init_train)], ret_forecast = stock_returns[i,] - risk_frees[-c(1:init_train)],
                                                volatility_forecast = vol_forecasts, gamma = 3, c_bp = 50)
}

for (i in c((length(stock_returns[,1]) +1):(2*length(stock_returns[,1])))){
  asset_management_sim[[i]] <- asset_allocation(actual_excess_rets = actual_stock_returns[-c(1:(init_train+init_valid))] - risk_frees[-c(1:(init_train+init_valid))],
                                                risk_frees = risk_frees[-c(1:(init_train+init_valid))], ret_forecast = stock_returns_val[i-length(stock_returns[,1]),] - risk_frees[-c(1:(init_train+init_valid))],
                                                volatility_forecast = vol_forecasts[-c(1:init_valid)], gamma = 3, c_bp = 50)
}

asset_management_sim[[length(asset_management_sim) + 1]] <- asset_allocation(actual_excess_rets = actual_stock_returns[-c(1:init_train)] - risk_frees[-c(1:init_train)],
                                                                             risk_frees = risk_frees[-c(1:init_train)], ret_forecast = hist_avg_forecasts - risk_frees[-c(1:init_train)],
                                                                             volatility_forecast = vol_forecasts, gamma = 3, c_bp = 50)

asset_management_sim[[length(asset_management_sim) + 1]] <- asset_allocation(actual_excess_rets = actual_stock_returns[-c(1:(init_train+init_valid))] - risk_frees[-c(1:(init_train+init_valid))],
                                                                             risk_frees = risk_frees[-c(1:(init_train+init_valid))], ret_forecast = hist_avg_forecasts[-c(1:init_valid)] - risk_frees[-c(1:(init_train+init_valid))],
                                                                             volatility_forecast = vol_forecasts[-c(1:init_valid)], gamma = 3, c_bp = 50)



# rename lists by predictive model name
names(asset_management_sim) <- c("trad","best_ncl","best_campbell","best_noneg","best_coef",
                                 "trad_val","ncl_val","campbell_val","nonneg_val","coef_val", "hist_avg", "hist_avg_val")


# calculate yearly utility gain over hist average

for (pred_model in c(1:length(asset_management_sim))){
  if (pred_model < 6 | pred_model == 11){
    asset_management_sim[[pred_model]]$yearly_util_gain <- (asset_management_sim[[pred_model]]$avg_utility - 
                                                              asset_management_sim$hist_avg$avg_utility)*400
    
  } else if (pred_model < 11 | pred_model == 12){
    asset_management_sim[[pred_model]]$yearly_util_gain <- (asset_management_sim[[pred_model]]$avg_utility - 
                                                              asset_management_sim$hist_avg_val$avg_utility)*400
  } 
}



### OOS EVALUATION


### calculate oos R^2-s

R_squared <- matrix(nrow=4, ncol = num_of_models)

for (j in c(1:nrow(R_squared))){
  for (i in c(1:num_of_models)){
    if (j == 1){
      R_squared[j,i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts,
                                       predictions = as.vector(univariate_models[[i]][[1]]))
    } else if (j == 2){
      R_squared[j,i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts,
                                       predictions = as.vector(univar_campbell[[i]][[1]]))
    } else if (j == 3){
      R_squared[j,i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts,
                                       predictions = as.vector(univar_nonneg[[i]][[1]]))
    } else if (j == 4){
      R_squared[j,i] <- calc_r_squared(targets= eq_prem_oos, hist_avg_pred = hist_avg_forecasts,
                                       predictions = as.vector(univar_coef[[i]][[1]]))
    }
    
  }
}

R_squared <- R_squared*100



## clark-west test

cw_p_vals <- vector(length = num_of_models)

for (i in c(1:num_of_models)){
  cw_p_vals[i] <- clark_west_test(targets = eq_prem_oos, predictions = univariate_models[[i]][[1]],
                                  hist_avg_pred = hist_avg_forecasts)[2]
  
}


plotting <- data.frame(lambdas = seq(from = 0.1, to = 1, by = 0.1),
                       R_squared_unresc = R_squared[1,-1], 
                       R_squared_campbell = R_squared[2,-1],
                       R_squared_nonneg = R_squared[3,-1],
                       R_squared_coef = R_squared[4,-1])


R_squared_plot1 <- ggplot(plotting, aes(x= lambdas, y = R_squared_unresc)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  geom_point(data=plotting[which.max(plotting$R_squared_unresc),], aes(x=lambdas, y=R_squared_unresc), colour="red", size=5) +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2") 

R_squared_plot2 <- ggplot(plotting, aes(x= lambdas, y = R_squared_campbell)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  geom_point(data=plotting[which.max(plotting$R_squared_campbell),], aes(x=lambdas, y=R_squared_campbell), colour="red", size=5) +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2") 

R_squared_plot3 <- ggplot(plotting, aes(x= lambdas, y = R_squared_nonneg)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  geom_point(data=plotting[which.max(plotting$R_squared_nonneg),], aes(x=lambdas, y=R_squared_nonneg), colour="red", size=5) +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2 (%)") 

R_squared_plot4 <- ggplot(plotting, aes(x= lambdas, y = R_squared_coef)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -2, to = 5.5, by = 0.5), limits = c(-2,5.5)) +
  geom_point(data=plotting[which.max(plotting$R_squared_coef),], aes(x=lambdas, y=R_squared_coef), colour="red", size=5) +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2 (%)") 



grid.arrange(R_squared_plot1, R_squared_plot2, R_squared_plot3, R_squared_plot4, ncol=2)


R_squared_plot1 <- ggplot(plotting, aes(x= lambdas, y = R_squared_unresc)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -2, to = 5.5, by = 0.5), limits = c(0,5.5)) +
  geom_point(data=plotting[which.max(plotting$R_squared_unresc),], aes(x=lambdas, y=R_squared_unresc), colour="red", size=5) +
  geom_hline(yintercept=R_squared[1,1], linetype="dashed") +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2 (%)") +
  ggtitle("Unrestricted") +
  theme(plot.title = element_text(hjust = 0.5))

R_squared_plot2 <- ggplot(plotting, aes(x= lambdas, y = R_squared_campbell)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -2, to = 5.5, by = 0.5), limits = c(0,5.5)) +
  geom_point(data=plotting[which.max(plotting$R_squared_campbell),], aes(x=lambdas, y=R_squared_campbell), colour="red", size=5) +
  geom_hline(yintercept=R_squared[2,1], linetype="dashed") +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2 (%)") +
  ggtitle("Nonneg & Coef Restrictions") +
  theme(plot.title = element_text(hjust = 0.5))

R_squared_plot3 <- ggplot(plotting, aes(x= lambdas, y = R_squared_nonneg)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -2, to = 5.5, by = 0.5), limits = c(0,5.5)) +
  geom_point(data=plotting[which.max(plotting$R_squared_nonneg),], aes(x=lambdas, y=R_squared_nonneg), colour="red", size=5) +
  geom_hline(yintercept=R_squared[3,1], linetype="dashed") +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2 (%)") +
  ggtitle("Nonneg Restriction") +
  theme(plot.title = element_text(hjust = 0.5)) 

R_squared_plot4 <- ggplot(plotting, aes(x= lambdas, y = R_squared_coef)) + 
  geom_line() +
  geom_point(color="blue", size = 3, shape = 17) +
  scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -2, to = 5.5, by = 0.5), limits = c(0,5.5)) +
  geom_point(data=plotting[which.max(plotting$R_squared_coef),], aes(x=lambdas, y=R_squared_coef), colour="red", size=5) +
  geom_hline(yintercept=R_squared[4,1], linetype="dashed") +
  xlab("Lambda values") +
  ylab(" Out of Sample R^2") +
  ggtitle("Coef Restriction") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(R_squared_plot1, R_squared_plot2, R_squared_plot3, R_squared_plot4, ncol=2)






hist_avg_forecasts <- forecast_hist_avg(y_data = eq_prem, init_train = init_train)
