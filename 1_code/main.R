# set working directory if necessary first
setwd("C:/Users/user/Desktop/git/tdk_2021")

source("1_code/project_functions.R")

# get libraries
library(tictoc)
library(abind)
library(ggplot2)
library(gridExtra)

# set storage list
X <- list()
X$DGPs <- list()


# set changing params
X$changing_params <- list()

r_squareds <- c(0.01,0.025,0.05,0.1,0.25)
correls <- c(0.0,0.4,0.8)

#r_squareds <- c(0.01)
#correls <- c(0.2)

X$changing_params$r_squareds <- r_squareds
X$changing_params$correls <- correls

len <- 300
sim_runs <- 100
i <- 1



# simulate DGPS with changing params
for (correl in X$changing_params$correls){
  for (r_squared in X$changing_params$r_squareds){
    X$DGPs <- def_dgp_params(X$DGPs, coefs = c(rep(1, times = 10)), r_squared = r_squared, correl = correl, 
                        num_of_vars = 10, len = len, sim_runs = sim_runs)
    X$DGPs <- simulate_dgp(X$DGPs,i, covar = "constant")
    
    i <- i + 1
  }
}


# train ncl methods
init_train <- 21
lambda <- seq(from = 0.1, to = 1.0, by = 0.1)

for (penalty_param in lambda){
  X$DGPs <- sim_forecasts(X$DGPs, init_train = init_train, method = "ncl", lambda = penalty_param)
}

# train OLS fc
X$DGPs <- sim_forecasts(X$DGPs, init_train = init_train, method = "traditional")



# estimate optimal weights

X$DGPs  <- ow_estim_dgps(X$DGPs, forecast_start = 71, weight_estim_ratio_vector = c(0.3,0.5,0.7), non_zero_weights = 0)

X$DGPs  <- ow_estim_dgps(X$DGPs, forecast_start = 71, weight_estim_ratio_vector = c(0.3,0.5,0.7), non_zero_weights = 1)



# get results

# dimensions in order: (R_squared, correl, coef, model_number, sim_run, data_ponint)

coefs <- X$DGPs$DGP_1$dgp_param$coefs

res_array <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                         coefs = coefs, method = "fc")

y_array <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                       coefs = coefs, method = "targets")

MSE_array <- get_MSE(res_array, y_array)


# normalize by OLS results


MSE_array_normalized <- MSE_normalize(MSE_array, OLS_index = 11)
MSE_array_normalized_nOLS <- MSE_array_normalized[,,,c(1:10),]


# MSE plots

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(apply(MSE_array_normalized_nOLS[r_squared,correl,,], 1, mean), 
         xaxt = "n", xlab ="NCL with lambda = 0.1, 0.2, ...",
         ylab = "NCL/OLS ratio", main = paste0("R^2 = ", as.character(r_squareds[r_squared]*100), 
                                               "%, cor param = ", as.character(correls[correl])))
    axis(1, at = seq(1, 10, by = 1))
    
  }
}

# MSE plots with kitchen sink

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(apply(MSE_array[r_squared,correl,,,], 1, mean), 
         xaxt = "n", xlab ="NCL with lambda = 0.1, 0.2, ...",
         ylab = "NCL/OLS ratio", main = paste0("R^2 = ", as.character(r_squareds[r_squared]*100), 
                                               "%, cor param = ", as.character(correls[correl])))
    axis(1, at = seq(1, 10, by = 1))
    
  }
}




















# get results for OWs

#res_array_unresc_OW <- get_results(X$DGPs, coefs = coefs, method = "OW", ow_type_index = 5)

MSE_array_unresc_OW <- get_MSE(res_array_unresc_OW, y_array, method = "OW")

res_array_nonneg_OW <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                                   coefs = coefs, method = "OW", ow_type_index = 5)

MSE_array_nonneg_OW <- get_MSE(res_array_nonneg_OW, y_array, method = "OW")



# Nonneq OW MSE plots

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(apply(MSE_array_nonneg_OW[r_squared, correl, , c(1:8,11), ,3], 1, mean), 
         xaxt = "n", xlab ="NCL with lambda = 0.1, 0.2, ...",
         ylab = "NCL/OLS ratio", main = paste0("R^2 = ", as.character(r_squareds[r_squared]), 
                                               "%, cor param = ", as.character(correls[correl])))
    
    
    
  }
}


# EW and OW MSE diffs

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(apply(MSE_array_nonneg_OW[r_squared, correl, , c(1:4,11), ,3], 1, mean) -  apply(MSE_array[r_squared,correl,,c(1:4,11),], 1, mean), 
         xaxt = "n", xlab ="NCL with lambda = 0.1, 0.2, ...",
         ylab = "NCL/OLS ratio", main = paste0("R^2 = ", as.character(r_squareds[r_squared]*100), 
                                               "%, cor param = ", as.character(correls[correl])))
    axis(1, at = seq(1, 22, by = 1))
    
  }
}


## magy változók egymáson való reg-jébol R^2 mindegyik idosorra

r_sq <- get_ow_r_squareds(X$DGPs, correls = correls, r_squareds = r_squareds,
                          coefs = coefs)

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(r_sq[r_squared, correl, 1,], 
         xaxt = "n", xlab ="NCL with lambda = 0.1, 0.2, ...",
         ylab = "R^2 of Var1 on Var2", main = paste0("R^2 = ", as.character(r_squareds[r_squared]*100), 
                                                     "%, cor param = ", as.character(correls[correl])))
    axis(1, at = seq(1, 22, by = 1))
    
  }
}

# check correlations

r_squareds <- c(0.01,0.025,0.05,0.1,0.25)
correls <- c(0.0,0.4,0.8)
correl_array <- correl_calc(X$DGPs, coefs = coefs)

# correl checks for weight estim regression

# cor(X1,X2)
par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,1])
  }
}

# cor(X1,X10)

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,2])
  }
}


# cor(X2,X10)

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,3])
  }
}




par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,1]+correl_array[r_squared, correl, 1,,2])
  }
}

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,4])
  }
}

# add up abs values of correlation

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,5])
  }
}

# correl values from "X1-XP" regressions

par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,6])
  }
}


par(mfrow=c(3,5))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(correl_array[r_squared, correl, 1,,7])
  }
}

### weights histogram

hist_array <- create_hist_array(X$DGPs, coefs = coefs, ow_index = 5, sim_num_max = 100, num_var_max = 1)


# weight hist for same model with diff params

par(mfrow=c(3,5))

model_num <- 10
for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    hist(hist_array[r_squared, correl, 1,model_num,1,,2])
    
  }
}

# wieght hist for lambda = 1 and OLS

par(mfrow=c(2,5))

model_nums <- c(10,11)
correl <- 3
for (model_num in model_nums){
  for (r_squared in c(1:5)){
    hist(hist_array[r_squared, correl, 1,model_num,1,,2])
    
  }
}


# weight hists

par(mfrow=c(3,11))

correl <- 2
for (r_squared in c(1,3,5)){
  for (model in c(1:11)){
    hist(hist_array[r_squared, correl, 1,model,1,,1])
  }
}


# standard deviation of weights

par(mfrow=c(3,5))

num_of_var <- 1
for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    plot(apply(hist_array[r_squared,correl,1,,num_of_var,,1],c(1),sd))
  }
}


sd_matrix <- matrix(data = NA, nrow = 5, ncol = 11)
correl <- 2
for (r_squared in c(1:5)){
  sd_matrix[r_squared, ] <- apply(hist_array[r_squared,correl,1,,num_of_var,,1],c(1),sd)
}
sd_matrix

cor_mat_array <- get_cor_mat(X$DGPs, coefs = coefs, method = "normal")
correl <- 2
model_num <- 11
a <- matrix(data = cor_mat_array[3,correl,1,model_num,], nrow = 10)
model_num <- 10
b <- matrix(data = cor_mat_array[3,correl,1,model_num,], nrow = 10)

cor_mat_array_min <- get_cor_mat(X$DGPs, coefs = coefs, method = "minus")
correl <- 2
model_num <- 11
a_2 <- matrix(data = cor_mat_array_min[3,correl,1,model_num,], nrow = 9)
model_num <- 10
b_2 <- matrix(data = cor_mat_array_min[3,correl,1,model_num,], nrow = 9)




# rm(list = ls())
# load("C:/Users/user/Desktop/git/tdk_2021/.RData")

#for (dgp_num in c(1:15)){
#  X[[1]][[dgp_num]][[3]] <- NULL
#}

# train p_ncl models

for (penalty_param in lambda){
  X$DGPs <- sim_forecasts(X$DGPs, init_train = init_train, method = "p_ncl", lambda = penalty_param)
}

# train kitchen sink
X$DGPs <- sim_forecasts(X$DGPs, init_train = init_train, method = "traditional", models = "kitchen_sink")


# validation of lambda

X$DGPs <- validate_lambda(X$DGPs, 71, 30, method = "expanding")

res_validated <- get_validated_res(X$DGPs)

y_array_valid <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                       coefs = coefs, method = "targets", fc_length = ncol(X$DGPs$DGP_1$validated_lambda[[1]]))

MSE_array_valid <- get_MSE(res_validated, y_array_valid)

#nonvalidated data

res_array_nv <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                         coefs = coefs, method = "fc", fc_length = ncol(X$DGPs$DGP_1$validated_lambda[[1]]))

y_array_nv <- get_results(X$DGPs, correls = correls, r_squareds = r_squareds,
                       coefs = coefs, method = "targets", 
                       fc_length = ncol(X$DGPs$DGP_1$validated_lambda[[1]]))

MSE_array_nv <- get_MSE(res_array_nv, y_array_nv)


# plot MSE

MSE_plots <- data.frame(r_squared_0.01_correl_0 = apply(MSE_array_nv[1,1,,c(1:10),], 1, mean))
MSE_ols_univar <- data.frame(r_squared_0.01_correl_0 = apply(MSE_array_nv[1,1,,c(1:11),], 1, mean)[11])
MSE_plots[[1]] <- MSE_plots[[1]] / MSE_ols_univar[[1]]

MSE_val_for_plots <- data.frame(r_squared_0.01_correl_0 = apply(MSE_array_valid[1,1,,,,drop=F], 1, mean))
MSE_val_for_plots[[1]] <-  MSE_val_for_plots[[1]] / MSE_ols_univar[[1]]

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    if (correl != 1 | r_squared != 1){
      MSE_plots[[length(MSE_plots) + 1]] <- apply(MSE_array_nv[r_squared,correl,,c(1:10),], 1, mean)
      MSE_ols_univar[[length(MSE_ols_univar) + 1]] <- apply(MSE_array_nv[r_squared,correl,,c(1:11),], 1, mean)[11]
      MSE_plots[[length(MSE_plots)]] <- MSE_plots[[length(MSE_plots)]] / MSE_ols_univar[[length(MSE_ols_univar)]]
      
      MSE_val_for_plots[[length(MSE_val_for_plots)+1]] <- apply(MSE_array_valid[r_squared,correl,,,,drop=F], 1, mean)
      MSE_val_for_plots[[length(MSE_val_for_plots)]] <- MSE_val_for_plots[[length(MSE_val_for_plots)]] / MSE_ols_univar[[length(MSE_ols_univar)]]
      
      }
  }
}




colnames(MSE_plots) <- c(colnames(MSE_plots)[1], "r_squared_0.025_correl_0", 
                         "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                         "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                         "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                         "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                         "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                         "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                         "r_squared_0.25_correl_0.8")
colnames(MSE_ols_univar) <- c(colnames(MSE_ols_univar)[1], "r_squared_0.025_correl_0", 
                         "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                         "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                         "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                         "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                         "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                         "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                         "r_squared_0.25_correl_0.8")

colnames(MSE_val_for_plots) <- c(colnames(MSE_val_for_plots)[1], "r_squared_0.025_correl_0", 
                              "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                              "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                              "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                              "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                              "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                              "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                              "r_squared_0.25_correl_0.8")

MSE_plots$numbers <- c(1:10)

MSE_plot_list <- list()

# create title vector
r_squareds <- c("1%", "2.5%", "5%", "10%", "25%")
correls <- c("0","0.4","0.8")
title_vector <- vector(length = 15)
i <- 1
for (correl in correls){
  for (r_squared in r_squareds){
    title_vector[i] <- paste0("R^2 = ", r_squared, ", Sigma = ", correl)
    i <- i + 1
    }
}


for (i in c(1:15)){
  MSE_plot_list[[length(MSE_plot_list)+1]] <- ggplot(MSE_plots, aes_string(x="numbers", y = names(MSE_plots)[i])) +
    geom_line() +
    geom_point(size = 3, shape = 17, colour = "blue") +
    scale_x_continuous(breaks = seq(from = 2, to = 10, by = 2), labels = seq(from = 0.2, to = 1, by = 0.2), limits = c(0,10)) +
    geom_point(x=which.min(MSE_plots[[i]]), y = MSE_plots[[i]][[which.min(MSE_plots[[i]])]], colour = "red", size = 4) +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_hline(yintercept=MSE_val_for_plots[[i]], linetype="solid", colour = "red") +
    xlab("Lambda values") +
    ylab("MSE ratio") +
    ggtitle(title_vector[i]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
}

grid.arrange(MSE_plot_list[[1]], MSE_plot_list[[2]], MSE_plot_list[[3]], MSE_plot_list[[4]],
             MSE_plot_list[[5]], MSE_plot_list[[6]], MSE_plot_list[[7]], MSE_plot_list[[8]],
             MSE_plot_list[[9]], MSE_plot_list[[10]], MSE_plot_list[[11]], MSE_plot_list[[12]],
             MSE_plot_list[[13]], MSE_plot_list[[14]], MSE_plot_list[[15]], ncol=5)


grid.arrange(MSE_plot_list[[6]], MSE_plot_list[[8]], MSE_plot_list[[10]], ncol = 3)




ggplot(MSE_plots, aes_string(x="numbers", y = names(MSE_plots)[i])) +
  geom_line() +
  geom_point(size = 3, shape = 17, colour = "blue") +
  scale_x_continuous(breaks = c(1:10), labels = seq(from = 0.1, to = 1, by = 0.1)) +
  geom_point(x=which.min(MSE_plots[[i]]), y = MSE_plots[[i]][[which.min(MSE_plots[[i]])]], colour = "red", size = 4) +
  geom_hline(yintercept=1, linetype="dashed") +
  xlab("Lambda values") +
  ylab("MSE") +
  ggtitle("R^2 = 1%, Sigma = 0") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(MSE_plots, aes(x=c(1:10), y = r_squared_0.025_correl_0)) +
  geom_line() +
  geom_point(size = 3, shape = 17, colour = "blue") +
  scale_x_continuous(breaks = c(1:10), labels = seq(from = 0.1, to = 1, by = 0.1)) +
  geom_point(x=which.min(MSE_plots$r_squared_0.025_correl_0), y = MSE_plots$r_squared_0.025_correl_0[[which.min(MSE_plots$r_squared_0.01_correl_0)]], colour = "red", size = 4) +
  geom_hline(yintercept=1, linetype="dashed") +
  xlab("Lambda values") +
  ylab("MSE") +
  ggtitle("R^2 = 1%, Sigma = 0") +
  theme(plot.title = element_text(hjust = 0.5))
  



#### Estimate bias-var-covar

comp_matrix <- matrix(data = NA, nrow = 11, ncol = 4)

predictors_list <- list()

for (i in c(1:11)){
  for (j in c(1:10)){
    
    if (j == 1){
      predictors_list[[i]] <- list()
    }
    predictors_list[[i]][[j]] <- X$DGPs$DGP_8$forecast_results[[i]][[1]][[3]][j,]
  }
}

for (i in c(1:11)){
  comp_matrix[i,1:3] <- calc_bias_var_covar(predictions_list = predictors_list[[i]], targets = X$DGPs$DGP_8$sim_data$y[c(22:300)])
}

comp_matrix[,4] <- rowSums(comp_matrix[,1:3])

### Calc lambda vals for different dgps


lambdas_for_params_val <- data.frame(r_squared_0.01_correl_0 = vector(length = 10))
i <- 1

for (dgp in X$DGPs){
  lambda_vector <- vector(length = 10)
  for (sim_run in c(1:100)){
    for (lambda in c(1:10)){
      lambda_vector[lambda] <- lambda_vector[lambda] + 
        sum(dgp$validated_lambda[[sim_run]][2,][which(dgp$validated_lambda[[sim_run]][2,] == lambda/10)])/(lambda/10)
    }
  }
  lambdas_for_params_val[[i]] <- lambda_vector
  i <- i +1
}

lambdas_for_params_val <- lambdas_for_params_val/sum(lambdas_for_params_val[[1]])

lambdas_for_params_val[[length(lambdas_for_params_val) +1]] <- c(1:10)

colnames(lambdas_for_params_val) <- c(colnames(lambdas_for_params_val)[1], "r_squared_0.025_correl_0", 
                                 "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                                 "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                                 "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                                 "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                                 "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                                 "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                                 "r_squared_0.25_correl_0.8", "numbers")

### LAMBDA PLOTS

val_lambdas_plot_list <- list()

for (i in c(1:15)){
  val_lambdas_plot_list[[i]] <- ggplot(data = lambdas_for_params_val, aes_string(x="numbers", y = names(lambdas_for_params_val)[i])) +
    geom_line() +
    geom_point(size = 3, shape = 17, colour = "blue") +
    geom_point(x=which.min(MSE_plots[[i]]), y = lambdas_for_params_val[[i]][[which.min(MSE_plots[[i]])]], colour = "red", size = 4) +
    scale_x_continuous(breaks = c(1:10), labels = seq(from = 0.1, to = 1, by = 0.1)) +
    scale_y_continuous(breaks = c(1:10)/10, labels = seq(from = 0.1, to = 1, by = 0.1), limits = c(0,1)) +
    xlab("Lambda values") +
    ylab("Relative frequency") +
    ggtitle(title_vector[1]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size =8))
}

grid.arrange(val_lambdas_plot_list[[1]], val_lambdas_plot_list[[2]], val_lambdas_plot_list[[3]], 
             val_lambdas_plot_list[[4]], val_lambdas_plot_list[[5]], val_lambdas_plot_list[[6]],
             val_lambdas_plot_list[[7]], val_lambdas_plot_list[[8]], val_lambdas_plot_list[[9]],
             val_lambdas_plot_list[[10]], val_lambdas_plot_list[[11]], val_lambdas_plot_list[[12]],
             val_lambdas_plot_list[[13]], val_lambdas_plot_list[[14]], val_lambdas_plot_list[[15]],
             ncol = 5)

ggplot(data = lambdas_for_params_val, aes(x=c(1:10), y = r_squared_0.01_correl_0)) +
  geom_line() +
  geom_point(size = 3, shape = 17, colour = "blue") +
  geom_point(x=which.min(MSE_plots[[1]]), y = lambdas_for_params_val[[1]][[which.min(MSE_plots[[1]])]], colour = "red", size = 4) +
  scale_x_continuous(breaks = c(1:10), labels = seq(from = 0.1, to = 1, by = 0.1)) +
  xlab("Lambda values") +
  ylab("Relative frequency") +
  ggtitle(title_vector[1]) +
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size =8))
  

### OW plots

# Nonneq OW MSE plots

OW_MSE_plots <- data.frame(r_squared_0.01_correl_0 = apply(MSE_array_nonneg_OW[1,1,,c(1:10),,2], 1, mean))
OW_ols_univar <- data.frame(r_squared_0.01_correl_0 = apply(MSE_array_nonneg_OW[1,1,,c(1:11), ,2], 1, mean)[11])
OW_MSE_plots[[1]] <- OW_MSE_plots[[1]] / OW_ols_univar[[1]]

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    if (correl != 1 | r_squared != 1){
      OW_MSE_plots[[length(OW_MSE_plots) + 1]] <- apply(MSE_array_nonneg_OW[r_squared,correl,,c(1:10),,2], 1, mean)
      OW_ols_univar[[length(OW_ols_univar) + 1]] <- apply(MSE_array_nonneg_OW[r_squared,correl,,c(1:11), ,2], 1, mean)[11]
      OW_MSE_plots[[length(OW_MSE_plots)]] <- OW_MSE_plots[[length(OW_MSE_plots)]] / OW_ols_univar[[length(OW_ols_univar)]]
    }
  }
}

OW_MSE_plots$numbers <- c(1:10)

colnames(OW_MSE_plots) <- c(colnames(OW_MSE_plots)[1], "r_squared_0.025_correl_0", 
                         "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                         "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                         "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                         "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                         "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                         "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                         "r_squared_0.25_correl_0.8", "numbers")
colnames(OW_ols_univar) <- c(colnames(OW_ols_univar)[1], "r_squared_0.025_correl_0", 
                              "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                              "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                              "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                              "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                              "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                              "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                              "r_squared_0.25_correl_0.8")

OW_MSE_plot_list <- list()


for (i in c(1:15)){
  OW_MSE_plot_list[[length(OW_MSE_plot_list)+1]] <- ggplot(OW_MSE_plots, aes_string(x="numbers", y = names(OW_MSE_plots)[i])) +
    geom_line() +
    geom_point(size = 3, shape = 17, colour = "blue") +
    scale_x_continuous(breaks = seq(from = 2, to = 10, by = 2), labels = seq(from = 0.2, to = 1, by = 0.2), limits = c(0,10)) +
    geom_point(x=which.min(OW_MSE_plots[[i]]), y = OW_MSE_plots[[i]][[which.min(OW_MSE_plots[[i]])]], colour = "red", size = 4) +
    geom_hline(yintercept=1, linetype="dashed") +
    xlab("Lambda values") +
    ylab("MSE ratio") +
    ggtitle(title_vector[i]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
}

grid.arrange(OW_MSE_plot_list[[1]], OW_MSE_plot_list[[2]], OW_MSE_plot_list[[3]], OW_MSE_plot_list[[4]],
             OW_MSE_plot_list[[5]], OW_MSE_plot_list[[6]], OW_MSE_plot_list[[7]], OW_MSE_plot_list[[8]],
             OW_MSE_plot_list[[9]], OW_MSE_plot_list[[10]], OW_MSE_plot_list[[11]], OW_MSE_plot_list[[12]],
             OW_MSE_plot_list[[13]], OW_MSE_plot_list[[14]], OW_MSE_plot_list[[15]], ncol=5)



### MSE diffs

MSE_diffs <- data.frame(r_squared_0.01_correl_0 = ((apply(MSE_array_nonneg_OW[1,1,,c(1:10),,2], 1, mean) - apply(MSE_array[1,1,,c(1:10),], 1, mean)) - 
                                                     (apply(MSE_array_nonneg_OW[1,1,,c(1:11),,2], 1, mean)[11] - apply(MSE_array[1,1,,c(1:11),], 1, mean)[11])))

for (correl in c(1:3)){
  for (r_squared in c(1:5)){
    if (correl != 1 | r_squared != 1){
      
      MSE_diffs[[length(MSE_diffs) + 1]] <- ((apply(MSE_array_nonneg_OW[r_squared,correl,,c(1:10),,2], 1, mean) - apply(MSE_array[r_squared,correl,,c(1:10),], 1, mean)) - 
                                               (apply(MSE_array_nonneg_OW[r_squared,correl,,c(1:11),,2], 1, mean)[11] - apply(MSE_array[r_squared,correl,,c(1:11),], 1, mean)[11]))
    }
    
  }
}

MSE_diffs <- MSE_diffs * 100

MSE_diffs$numbers <- c(1:10)

colnames(MSE_diffs) <- c(colnames(MSE_diffs)[1], "r_squared_0.025_correl_0", 
                            "r_squared_0.05_correl_0", "r_squared_0.1_correl_0",
                            "r_squared_0.25_correl_0", "r_squared_0.01_correl_0.4",
                            "r_squared_0.025_correl_0.4", "r_squared_0.05_correl_0.4",
                            "r_squared_0.1_correl_0.4", "r_squared_0.25_correl_0.4",
                            "r_squared_0.01_correl_0.8", "r_squared_0.025_correl_0.8",
                            "r_squared_0.05_correl_0.8", "r_squared_0.1_correl_0.8",
                            "r_squared_0.25_correl_0.8", "numbers")




MSE_diffs_plot_list <- list()


for (i in c(1:15)){
  MSE_diffs_plot_list[[length(MSE_diffs_plot_list)+1]] <- ggplot(MSE_diffs, aes_string(x="numbers", y = names(MSE_diffs)[i])) +
    geom_line() +
    geom_point(size = 3, shape = 17, colour = "blue") +
    scale_x_continuous(breaks = seq(from = 2, to = 10, by = 2), labels = seq(from = 0.2, to = 1, by = 0.2), limits = c(0,10)) +
    geom_point(x=which.min(MSE_diffs[[i]]), y = MSE_diffs[[i]][[which.min(MSE_diffs[[i]])]], colour = "red", size = 4) +
    xlab("Lambda values") +
    ylab("MSE diff") +
    ggtitle(title_vector[i]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
}

grid.arrange(MSE_diffs_plot_list[[1]], MSE_diffs_plot_list[[2]], MSE_diffs_plot_list[[3]], MSE_diffs_plot_list[[4]],
             MSE_diffs_plot_list[[5]], MSE_diffs_plot_list[[6]], MSE_diffs_plot_list[[7]], MSE_diffs_plot_list[[8]],
             MSE_diffs_plot_list[[9]], MSE_diffs_plot_list[[10]], MSE_diffs_plot_list[[11]], MSE_diffs_plot_list[[12]],
             MSE_diffs_plot_list[[13]], MSE_diffs_plot_list[[14]], MSE_diffs_plot_list[[15]], ncol=5)




