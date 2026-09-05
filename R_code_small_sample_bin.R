
rm(list = ls())

#Performance of G_computation estimators in high dimensional settings
library(mvtnorm)
library(dplyr)
library(lme4)
library(sandwich)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(gee)
library(stdReg)
library(drgee)
library(geepack)
library(nlme)
library(glmtoolbox)
library(metafor)
library(cowplot)
library(grid)
library(ggpubr)
library(speff2trial) 
library(readstata13)
library(LaplacesDemon)
library(xtable)
library(sandwich)
library(glmnet)
library(caret)
library(data.table)
library(arm)
library(MASS)
library(HOIFCar)
library(logistf)
##########################################################################################
#Functions
#Unadjusted 
Unadjusted <- function(data,trt, y, prob) {
  n <- dim(data)[1]
  tau_unadj <- mean(y[trt == 1]) - mean(y[trt == 0])
  se_var_unadj <- sqrt(1 / n * (var(y[trt==1]) / mean(trt) + var(y[trt==0]) / (1 - mean(trt))))
  return(list(tau_unadj, se_var_unadj))
}


#Standard G-computation estimators
stand_G_comp <- function(data,X,trt,y,prob) {
  cov <- c(names(data %>% dplyr::select(-y)))
  formula <- as.formula(paste("y~", paste(cov, collapse = "+")))
  model_fit <- bayesglm(formula, data = data,family = binomial(link = "logit"),x = TRUE)
  coef <- coef(model_fit)[2]
  Var_coef_HC1 <- vcovHC(model_fit, type = "HC1")[2, 2]
  SE_coef_HC1 <- sqrt(Var_coef_HC1)
  new_dta_0 <- data.frame(trt=0,X)
  new_dta_1 <- data.frame(trt=1,X)
  pred_0 <- as.numeric(predict(model_fit, newdata = new_dta_0,type="response"))
  pred_1 <- as.numeric(predict(model_fit, newdata = new_dta_1,type="response"))
  
  #TMLE step
  new_data <- data.frame(trt,X)
  pred <- as.numeric(predict(model_fit, newdata = new_data,type="response"))
  epsilon <- coef(glm(y~trt+offset(qlogis(pred)), family=binomial))
  Q01<-plogis(qlogis(pred_0)+epsilon[[1]]) 
  Q11<-plogis(qlogis(pred_1)+epsilon[[1]]+epsilon[[2]])
  G_comp <- mean(Q11-Q01)
  
  trt <- as.numeric(trt)
  inf <- as.vector(trt/prob*(y-Q11)+Q11-((1-trt)/(1-prob)*(y-Q01)+Q01))
  n <- dim(data)[1]
  p <-  dim(data)[2]-2
  se_G_comp_inf <- sqrt(var(inf)/n)
  correction_factor <- (n-1)/(n-p-2)
  se_G_comp_inf_corr <- sqrt(var(inf)/n*correction_factor)
  return(list(G_comp, se_G_comp_inf,se_G_comp_inf_corr,coef,SE_coef_HC1))
}

#G-computation estimators with LASSO
LASSO_G_comp <- function(data,X,trt,y,prob) {
  #Fit LASSO model and select covariates 
  X_main <- model.matrix(~ trt+X)[, -1] 
  penalty_factor <- ifelse(colnames(X_main) == "trt", 0, 1)
  lasso_model <- cv.glmnet(X_main, y, alpha = 1, family="binomial", penalty.factor = penalty_factor)
  LASSO_coef <- as.matrix(coef(lasso_model, s = "lambda.min"))
  non_zero_coef_vars <- rownames(LASSO_coef)[LASSO_coef!= 0]
  selected_vars <- non_zero_coef_vars[non_zero_coef_vars %in% colnames(X_main)]
  p_selected <- length(selected_vars)-1
  #Post-LASSO Regression and predictions
  X_selected <- X_main[, selected_vars, drop = FALSE]
  #Combine with outcome for post-LASSO regression
  data_selected <- data.frame(y, X_selected)
  if (!"trt" %in% names(data_selected)) {
    message("'trt' column not found.")
    data_selected$trt <- trt
  }
  cov <- names(data_selected %>% dplyr::select(-y))
  formula <- as.formula(paste("y ~", paste(cov, collapse = " + ")))
  model_fit <- bayesglm(formula, data = data_selected,family = binomial(link = "logit"),x = TRUE)
  coef <- coef(model_fit)[2]
  Var_coef_HC1 <- vcovHC(model_fit, type = "HC1")[2, 2]
  SE_coef_HC1 <- sqrt(Var_coef_HC1)
  new_dta_0 <- data.frame(X_selected)
  new_dta_1 <- data.frame(X_selected)
  new_dta_0[, "trt"] <- 0
  new_dta_1[, "trt"] <- 1
  pred_0 <- predict(model_fit, newdata = new_dta_0,type="response")
  pred_1 <- predict(model_fit, newdata = new_dta_1,type="response")
  
  #TMLE step
  new_data <- data.frame(X_selected)
  new_data[, "trt"] <- trt
  pred <- as.numeric(predict(model_fit, newdata = new_data,type="response"))
  epsilon <- coef(glm(y~trt+offset(qlogis(pred)), family=binomial))
  Q01<-plogis(qlogis(pred_0)+epsilon[[1]]) 
  Q11<-plogis(qlogis(pred_1)+epsilon[[1]]+epsilon[[2]])
  G_comp <- mean(Q11-Q01)
  
  inf <- as.vector(trt/prob*(y-Q11)+Q11-((1-trt)/(1-prob)*(y-Q01)+Q01)) 
  n <- dim(data)[1]
  se_G_comp_inf <- sqrt(var(inf)/n)
  #correction_factor
  correction_factor <- (n-1)/(n-p_selected-2)
  se_G_comp_inf_corr <- sqrt(var(inf)/n*correction_factor)
  return(list(G_comp,se_G_comp_inf,se_G_comp_inf_corr,p_selected,coef,SE_coef_HC1))
}

#G-computation estimators with cross fitting
CF_G_comp <- function(data,X,trt,y,prob,folds) {
  Q_0_all <- numeric(0)
  Q_1_all <- numeric(0)
  y_all <- numeric(0)
  trt_all <- numeric(0)
  for (i in seq_along(folds)) {
    test_indices <- folds[[i]]
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    #Train the model on training data set
    cov <- c(names(train_data %>% dplyr::select(-y)))
    formula <- as.formula(paste("y ~", paste(cov, collapse = " + ")))
    model_fit <- bayesglm(formula, data = train_data,family = binomial(link = "logit"),x = TRUE)
    
    #Predict on the test data
    test_data_0 <- test_data %>%dplyr::select(-y)
    test_data_0$trt <- 0  # Set treatment to 0
    test_data_1 <- test_data %>%dplyr::select(-y)
    test_data_1$trt <- 1  # Set treatment to 1
    
    #Predictions
    pred_0_i <- predict(model_fit, newdata = test_data_0,type="response")
    pred_1_i <- predict(model_fit, newdata = test_data_1,type="response")
    
    #TMLE step
    new_data <- test_data %>%dplyr::select(-y)
    pred_i <- as.numeric(predict(model_fit, newdata = new_data,type="response"))
    epsilon_i <- coef(glm(test_data$y~test_data$trt+offset(qlogis(pred_i)), family=binomial))
    Q01_i<-plogis(qlogis(pred_0_i)+epsilon_i[[1]]) 
    Q11_i<-plogis(qlogis(pred_1_i)+epsilon_i[[1]]+epsilon_i[[2]])
    
    Q_0_all <- c(Q_0_all, Q01_i)
    Q_1_all <- c(Q_1_all, Q11_i)
    y_all <- c(y_all, test_data$y)  
    trt_all <- c(trt_all, test_data$trt)  
  }
  
  G_comp <- mean(Q_1_all-Q_0_all)
  inf <- as.vector(trt_all/prob*(y_all-Q_1_all)+Q_1_all-((1-trt_all)/(1-prob)*(y_all-Q_0_all)+Q_0_all)) 
  n <- dim(data)[1]
  p <-  dim(data)[2]-2
  se_G_comp_inf <- sqrt(var(inf)/n)
  correction_factor <- (n-1)/(n-p-2)
  se_G_comp_inf_corr <- sqrt(var(inf)/n*correction_factor)
  return(list(G_comp,se_G_comp_inf,se_G_comp_inf_corr))
}

#G-computation estimators with cross fitting and LASSO
CF_LASSO_G_comp <- function(data,X,trt,y,prob,folds) {
  Q_0_all <- numeric(0)
  Q_1_all <- numeric(0)
  y_all <- numeric(0)
  trt_all <- numeric(0)
  for (i in seq_along(folds)) {
    test_indices <- folds[[i]]
    train_indices <- setdiff(seq_len(nrow(data)), test_indices)
    #Split train and test sets
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    #Ensure X, y, and trt are properly subsetted
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    trt_train <- trt[train_indices]
    # Ensure valid feature names
    X_train_df <- data.frame(trt_train, X_train)
    X_main <- model.matrix(~ ., data = X_train_df)[, -1]  
    colnames(X_main) <- make.names(colnames(X_main))
    penalty_factor <- ifelse(colnames(X_main) == "trt_train", 0, 1)
    # Fit LASSO model
    lasso_model <- cv.glmnet(X_main, y_train, alpha = 1, family="binomial", penalty.factor = penalty_factor)
    #Extract selected variables
    LASSO_coef <- as.matrix(coef(lasso_model, s = "lambda.min"))
    selected_vars <- rownames(LASSO_coef)[LASSO_coef != 0]
    selected_vars <- selected_vars[selected_vars %in% colnames(X_main)]
    #Ensure selected_vars contains at least `trt_train`
    if (length(selected_vars) == 0) {
      warning("No variables selected by LASSO. Using only trt.")
      selected_vars <- "trt_train"
    }
    #Fit regression
    X_selected <- X_main[, selected_vars, drop = FALSE]
    dta_selected <- data.frame(y_train, X_selected)
    if (!"trt_train" %in% names(dta_selected)) {
      dta_selected$trt_train <- trt_train
    }
    formula <- as.formula(paste("y_train ~", paste(names(dta_selected)[-1], collapse = " + ")))
    model_fit <- bayesglm(formula, data = dta_selected,family = binomial(link = "logit"),x = TRUE)
    # Handle missing or empty test data
    existing_vars <- selected_vars[selected_vars %in% colnames(test_data)]
    if (length(existing_vars) == 0) {
      warning("None of the selected variables exist in test_data. Using only trt.")
      test_data_0 <- data.frame(trt_train = rep(0, nrow(test_data)))
      test_data_1 <- data.frame(trt_train = rep(1, nrow(test_data)))
      new_data <- data.frame(trt_train = test_data$trt)
    } else {
      test_data_0 <- test_data[, existing_vars, drop = FALSE]
      test_data_1 <- test_data[, existing_vars, drop = FALSE]
      test_data_0$trt_train <- 0
      test_data_1$trt_train <- 1
      new_data <- test_data[, existing_vars, drop = FALSE]
      new_data$trt_train <- test_data$trt
    }
    colnames(test_data_0) <- make.names(colnames(test_data_0))
    colnames(test_data_1) <- make.names(colnames(test_data_1))
    colnames(new_data) <- make.names(colnames(new_data))
    
    #Predictions
    pred_0_i <- predict(model_fit, newdata = test_data_0,type="response")
    pred_1_i <- predict(model_fit, newdata = test_data_1,type="response")
    
    #TMLE step
    pred_i <- as.numeric(predict(model_fit, newdata = new_data,type="response"))
    epsilon_i <- coef(glm(test_data$y~test_data$trt+offset(qlogis(pred_i)), family=binomial))
    Q01_i<-plogis(qlogis(pred_0_i)+epsilon_i[[1]]) 
    Q11_i<-plogis(qlogis(pred_1_i)+epsilon_i[[1]]+epsilon_i[[2]])
    
    #Append to final vectors
    Q_0_all <- c(Q_0_all, Q01_i)
    Q_1_all <- c(Q_1_all, Q11_i)
    y_all <- c(y_all, y[test_indices])
    trt_all <- c(trt_all, trt[test_indices])
  }
  #Ensure all output vectors have the same length
  min_length <- min(length(y_all), length(trt_all), length(Q_0_all), length(Q_1_all))
  y_all <- y_all[seq_len(min_length)]
  trt_all <- trt_all[seq_len(min_length)]
  Q_0_all <- Q_0_all[seq_len(min_length)]
  Q_1_all <- Q_1_all[seq_len(min_length)]
  G_comp <- mean(Q_1_all-Q_0_all)
  inf <- as.vector(trt_all/prob*(y_all-Q_1_all)+Q_1_all-((1-trt_all)/(1-prob)*(y_all-Q_0_all)+Q_0_all))
  n <- dim(data)[1]
  se_G_comp_inf <- sqrt(var(inf)/n)
  return(list(G_comp,se_G_comp_inf))
}



#HOIF_motivated covariate adjusted estimators
HOIF_G_comp <- function(y,X,trt,prob){
  n <- nrow(X)
  
  ###tau_HOIF_1
  xc <- cbind(1,scale(X,scale=FALSE))
  H <- xc %*% ginv(t(xc) %*% xc) %*% t(xc)
  yt_adj2_hat <- as.numeric((H - diag(diag(H))) %*% (trt * y / prob))
  yc_adj2_hat <- as.numeric((H - diag(diag(H))) %*% ((1 - trt) * y / (1 - prob)))
  psi1_adj2_vec <- trt * y / prob-(trt/prob-1) * yt_adj2_hat
  psi0_adj2_vec <- (1 - trt) * y / (1 - prob)-((1 - trt)/(1 - prob)-1) * yc_adj2_hat
  
  tau1_adj2 <- mean(psi1_adj2_vec)
  tau0_adj2 <- mean(psi0_adj2_vec)
  tau_HOIF_1 <-tau1_adj2-tau0_adj2
  
  infl1_adj2 <- psi1_adj2_vec - tau1_adj2
  infl0_adj2 <- psi0_adj2_vec - tau0_adj2
  var_tau_HOIF_1 <- 1/n*var(infl1_adj2-infl0_adj2)
  
  
  ###HOIF_2
  tau1_unadj <- mean(diag(H)[trt==1]*y[trt==1])/mean(diag(H)[trt==1])
  tau0_unadj <- mean(diag(H)[trt==0]*y[trt==0])/mean(diag(H)[trt==0])
  
  yt_adj2c_hat <- as.numeric((H - diag(diag(H))) %*% (trt * (y - tau1_unadj) / prob))
  yc_adj2c_hat <- as.numeric((H - diag(diag(H))) %*% ((1 - trt) * (y - tau0_unadj) / (1 - prob)))
  
  psi1_adj2c_vec <- trt * y / prob-(trt/prob-1) * yt_adj2c_hat
  psi0_adj2c_vec <- (1 - trt) * y/(1 - prob)-((1 - trt) / (1 - prob)-1) * yc_adj2c_hat
  
  tau1_adj2c <- mean(psi1_adj2c_vec)
  tau0_adj2c <- mean(psi0_adj2c_vec)
  tau_HOIF_2 <- tau1_adj2c-tau0_adj2c
  
  
  infl1_adj2c <- trt * (y - tau1_adj2c) / prob + (1 - trt / prob) * yt_adj2c_hat
  infl0_adj2c <- (1 - trt) * (y - tau0_adj2c) / (1 - prob) + (1 - (1 - trt) / (1 - prob)) * yc_adj2c_hat
  
  var_tau_HOIF_2 <- 1/n*var(infl1_adj2c-infl0_adj2c)
  
  tau_vec <- c(tau_HOIF_1, tau_HOIF_2)
  names(tau_vec) <- c('HOIF_1', 'HOIF_2')
  var_infl_vec <- c(var_tau_HOIF_1, var_tau_HOIF_2)
  names(var_infl_vec) <- c('HOIF_1','HOIF_2')
  
  return(list(
    tau_vec = tau_vec,
    var_infl_vec = var_infl_vec))
  
}



# JASA estimator
fit_jasa <- function(y, X, trt, family_type) {
  
  # Center covariates without scaling by standard deviation
  Xc <- scale(X, scale = FALSE)
  
  # Add intercept column
  Xc_aug <- cbind(1, Xc)
  
  # Fit JASA estimator
  fit <- HOIFCar::fit.JASA(y, Xc_aug, trt, family = family_type, opt_obj = "beta")
  
  # Return uncalibrated and calibrated estimates with standard errors
  list(
    JASA_est     = fit[[1]][[1]],
    JASA_cal_est = fit[[1]][[2]],
    JASA_SE      = sqrt(fit[[4]][[1]]),
    JASA_cal_SE  = sqrt(fit[[4]][[2]])
  )
}

#Firth-based G-computation estimator for binary outcomes
firth_G_comp <- function(data, trt, X, y, prob) {
  
  # Build formula using all predictors except the outcome
  covars <- names(data)[names(data) != "y"]
  formula <- as.formula(paste("y ~", paste(covars, collapse = " + ")))
  
  # Fit Firth penalized logistic regression
  model_fit <- logistf::logistf(formula, data = data)
  
  # Counterfactual data sets under treatment = 0 and treatment = 1
  new_dta_0 <- data.frame(trt = 0, X)
  names(new_dta_0) <- names(data)[-1]
  
  new_dta_1 <- data.frame(trt = 1, X)
  names(new_dta_1) <- names(data)[-1]
  
  # Predict counterfactual probabilities
  pred_0 <- as.numeric(predict(model_fit, newdata = new_dta_0, type = "response"))
  pred_1 <- as.numeric(predict(model_fit, newdata = new_dta_1, type = "response"))
  
  # Predict observed-outcome probabilities under observed treatment
  new_data <- data.frame(trt, X)
  names(new_data) <- names(data)[-1]
  pred <- as.numeric(predict(model_fit, newdata = new_data, type = "response"))
  
  # Targeting step with offset logistic regression
  epsilon <- coef(glm(y ~ trt + offset(qlogis(pred)), family = binomial))
  
  # Updated predictions under treatment and control
  Q01 <- plogis(qlogis(pred_0) + epsilon[[1]])
  Q11 <- plogis(qlogis(pred_1) + epsilon[[1]] + epsilon[[2]])
  
  # Marginal treatment effect estimate
  tau_firth <- mean(pred_1 - pred_0)
  tau <- mean(Q11 - Q01)
  
  # Dimensions
  n_obs <- nrow(data)
  p_dim <- ncol(X)
  
  # Influence-function-based pseudo-outcome
  inf_firth<- as.vector(trt/prob * (y - pred_1) + pred_1 -((1 - trt) / (1 - prob) * (y - pred_0) + pred_0))
  inf <- as.vector(trt/prob * (y - Q11) + Q11 -((1 - trt) / (1 - prob) * (y - Q01) + Q01))
  
  # Standard error based on influence function
  se_firth <- sqrt(var(inf_firth)/n_obs)
  se_tau<- sqrt(var(inf)/n_obs)
  
  list(tau_firth = tau_firth,tau = tau,se_firth= se_firth, se_tau= se_tau)
}




Simulation_function <- function(n,q){
  #q <- 35
  #n <- 50
  
  #Generate X
  generate_cov_matrix <- function(q) {
    outer(1:q, 1:q, function(i, j) 0.1^abs(i-j))
  }
  
  Sigma <- generate_cov_matrix(q)
  mu <- rep(0, q)
  X <- mvrnorm(n, mu = mu, Sigma = Sigma) 
  
  #coefficients
  gamma <- (-1)^(1:q)/sqrt(1:q)
  
  #Linear predictor
  f_X <- X %*% gamma
  
#Treatment assignment
  trt <- rbinom(n,1,0.5) 
#Outcome
  pi <- invlogit(0.25+0.75*trt+f_X)
  y <- rbinom(n, 1, pi)
  data <- data.frame(y,trt,X)
  prob <- mean(trt)
  
  
##Unadjusted estimators
  unadjusted <- Unadjusted(data,trt,y,prob)
  tau_unadj <- unadjusted[[1]]
  se_tau_unadj_inf <- unadjusted[[2]]
  
  
##Standard G-computation estimators
  G_comp_stand <- stand_G_comp(data,X,trt,y,prob)
  tau_stand <- G_comp_stand[[1]]
  se_tau_stand_inf <- G_comp_stand[[2]]
  se_tau_stand_inf_corr <- G_comp_stand[[3]]
  coef_stand <- G_comp_stand[[4]]
  se_coef_stand <- G_comp_stand[[5]]
  
  
#G-computation estimators with LASSO
  G_comp_LASSO <- LASSO_G_comp(data,X,trt,y,prob)
  tau_LASSO <- G_comp_LASSO[[1]]
  se_tau_LASSO_inf <- G_comp_LASSO[[2]]
  se_tau_LASSO_inf_corr <- G_comp_LASSO[[3]]
  p_selected <- G_comp_LASSO[[4]]
  coef_LASSO <- G_comp_LASSO[[5]]
  se_coef_LASSO <- G_comp_LASSO[[6]]
  
#G-computation estimators with cross fitting
  fold_ids <- sample(rep(1:5, length.out = dim(data)[1]))
  folds <- split(1:dim(data)[1], fold_ids)
  G_comp_CF <- CF_G_comp(data,X,trt,y,prob,folds)
  tau_CF <- G_comp_CF[[1]]
  se_tau_CF_inf <- G_comp_CF[[2]]
  se_tau_CF_inf_corr <- G_comp_CF[[3]]
  
#G-computation estimators with cross fitting and LASSO
  G_comp_CF_LASSO <- CF_LASSO_G_comp(data,X,trt,y,prob,folds)
  tau_CF_LASSO <- G_comp_CF_LASSO[[1]]
  se_tau_CF_LASSO_inf <- G_comp_CF_LASSO[[2]]
  
  
#HOIF_motivated covariate adjusted estimators
  G_comp_HOIF <- HOIF_G_comp(y,X,trt,prob)
  tau_HOIF_1 <- G_comp_HOIF$tau_vec[[1]]
  se_tau_HOIF_1 <- sqrt(G_comp_HOIF$var_infl_vec[[1]])
  tau_HOIF_2 <- G_comp_HOIF$tau_vec[[2]]
  se_tau_HOIF_2 <- sqrt(G_comp_HOIF$var_infl_vec[[2]])
  
  
##JASA estimator
  JASA_estimator <- fit_jasa(y,X,trt,family_type="binomial")
  tau_JASA <- JASA_estimator$JASA_est
  se_tau_JASA <- JASA_estimator$JASA_SE
  tau_JASA_cal <- JASA_estimator$JASA_cal_est
  se_tau_JASA_cal <- JASA_estimator$JASA_cal_SE
  
  
##Firth-based G-computation estimator
  Firth_estimator <- firth_G_comp(data,trt,X,y,prob)
  tau_Firth <- Firth_estimator$tau_firth
  se_tau_Firth <- Firth_estimator$se_firth
  tau_Firth_target <- Firth_estimator$tau
  se_tau_Firth_target <- Firth_estimator$se_tau
  
  
  return(list(
    "tau_unadj"=tau_unadj,
    "se_tau_unadj_inf"=se_tau_unadj_inf,
    
    
    "tau_stand"=tau_stand,
    "se_tau_stand_inf"=se_tau_stand_inf,
    "se_tau_stand_inf_corr"=se_tau_stand_inf_corr,
    
    "tau_LASSO"=tau_LASSO,
    "se_tau_LASSO_inf"=se_tau_LASSO_inf,
    "se_tau_LASSO_inf_corr"=se_tau_LASSO_inf_corr,
    "p_selected"=p_selected,
    
    "tau_CF"=tau_CF,
    "se_tau_CF_inf"=se_tau_CF_inf,
    "se_tau_CF_inf_corr"=se_tau_CF_inf_corr,
    
    "tau_CF_LASSO"=tau_CF_LASSO,
    "se_tau_CF_LASSO_inf"=se_tau_CF_LASSO_inf,
    
    "tau_HOIF_1"=tau_HOIF_1,
    "se_tau_HOIF_1"=se_tau_HOIF_1,
    
    "tau_HOIF_2"=tau_HOIF_2,
    "se_tau_HOIF_2"=se_tau_HOIF_2,
    
    "tau_JASA"=tau_JASA,
    "se_tau_JASA"=se_tau_JASA,
    "tau_JASA_cal"=tau_JASA_cal,
    "se_tau_JASA_cal"=se_tau_JASA_cal,
    
    "tau_Firth"=tau_Firth,
    "se_tau_Firth"=se_tau_Firth,
    "tau_Firth_target"=tau_Firth_target,
    "se_tau_Firth_target"=se_tau_Firth_target))
}


#Specify simulation settings
simulate_settings <- function(nrep, setting_params) {
  results <- list()
  for (i in 1:length(setting_params)) {
    params <- setting_params[[i]]
    setting_results <- vector("list",length=nrep)
    for (j in 1:nrep) {
      setting_results[[j]] <- Simulation_function(params$n, params$q)
    }
    results[[i]] <- setting_results
  }
  return(results)
}

#q and n
#Settings
setting_params <- list(list(q=3, n=50),list(q=5, n=50),list(q=10, n=50),
                       list(q=15, n=50),list(q=20, n=50),list(q=35, n=50))


###Perform simulations
nrep <- 10
simulation_results <- simulate_settings(nrep, setting_params)

#Manage the simulated data sets
setting_data <- lapply(simulation_results, function(sim_list) {
  do.call(rbind, lapply(sim_list, as.data.frame))
})


#saveRDS(setting_data, "C:/Users/malene/OneDrive - UGent/Desktop/R_code_project_4/Simulated_data/bin/v4/Setting_1/bin_s1_sample_50.rds")
#setting_data <- readRDS("C:/Users/malene/OneDrive - UGent/Desktop/R_code_project_4/Simulated_data/bin/v4/Setting_1/bin_s1_sample_50.rds")

#Datasets
result_k_0.05_NA_Inf <- setting_data[[1]]
result_k_0.05_NA_Inf[result_k_0.05_NA_Inf == Inf|result_k_0.05_NA_Inf == -Inf] <- NA
result_k_0.05 <- na.omit(result_k_0.05_NA_Inf)

result_k_0.1_NA_Inf <- setting_data[[2]]
result_k_0.1_NA_Inf[result_k_0.1_NA_Inf == Inf|result_k_0.1_NA_Inf == -Inf] <- NA
result_k_0.1 <- na.omit(result_k_0.1_NA_Inf)

result_k_0.2_NA_Inf <- setting_data[[3]]
result_k_0.2_NA_Inf[result_k_0.2_NA_Inf == Inf|result_k_0.2_NA_Inf == -Inf] <- NA
result_k_0.2 <- na.omit(result_k_0.2_NA_Inf)

result_k_0.3_NA_Inf <- setting_data[[4]]
result_k_0.3_NA_Inf[result_k_0.3_NA_Inf == Inf|result_k_0.3_NA_Inf == -Inf] <- NA
result_k_0.3 <- na.omit(result_k_0.3_NA_Inf)

result_k_0.4_NA_Inf <- setting_data[[5]]
result_k_0.4_NA_Inf[result_k_0.4_NA_Inf == Inf|result_k_0.4_NA_Inf == -Inf] <- NA
result_k_0.4 <- na.omit(result_k_0.4_NA_Inf)

result_k_0.7_NA_Inf <- setting_data[[6]]
result_k_0.7_NA_Inf[result_k_0.7_NA_Inf == Inf|result_k_0.7_NA_Inf == -Inf] <- NA
result_k_0.7 <- na.omit(result_k_0.7_NA_Inf)
#####################################################################################

###True ATE
TE_k_0.05=0.1350419;TE_k_0.1=0.1299136;TE_k_0.2=0.1237476;TE_k_0.3=0.1204297;TE_k_0.4=0.1181031;TE_k_0.7=0.1140678


###Bias
###Unadjusted estimators
bias_tau_unadj_k_0.05 <- mean(result_k_0.05$tau_unadj-TE_k_0.05)
bias_tau_unadj_k_0.1 <- mean(result_k_0.1$tau_unadj-TE_k_0.1)
bias_tau_unadj_k_0.2 <- mean(result_k_0.2$tau_unadj-TE_k_0.2)
bias_tau_unadj_k_0.3 <- mean(result_k_0.3$tau_unadj-TE_k_0.3)
bias_tau_unadj_k_0.4 <- mean(result_k_0.4$tau_unadj-TE_k_0.4)
bias_tau_unadj_k_0.7 <- mean(result_k_0.7$tau_unadj-TE_k_0.7)

###Standard G-computation estimators
bias_tau_stand_k_0.05 <- mean(result_k_0.05$tau_stand-TE_k_0.05)
bias_tau_stand_k_0.1 <- mean(result_k_0.1$tau_stand-TE_k_0.1)
bias_tau_stand_k_0.2 <- mean(result_k_0.2$tau_stand-TE_k_0.2)
bias_tau_stand_k_0.3 <- mean(result_k_0.3$tau_stand-TE_k_0.3)
bias_tau_stand_k_0.4 <- mean(result_k_0.4$tau_stand-TE_k_0.4)
bias_tau_stand_k_0.7 <- mean(result_k_0.7$tau_stand-TE_k_0.7)

###G-computation estimators with LASSO
bias_tau_LASSO_k_0.05 <- mean(result_k_0.05$tau_LASSO-TE_k_0.05)
bias_tau_LASSO_k_0.1 <- mean(result_k_0.1$tau_LASSO-TE_k_0.1)
bias_tau_LASSO_k_0.2 <- mean(result_k_0.2$tau_LASSO-TE_k_0.2)
bias_tau_LASSO_k_0.3 <- mean(result_k_0.3$tau_LASSO-TE_k_0.3)
bias_tau_LASSO_k_0.4 <- mean(result_k_0.4$tau_LASSO-TE_k_0.4)
bias_tau_LASSO_k_0.7 <- mean(result_k_0.7$tau_LASSO-TE_k_0.7)

###G-computation estimators with CF
bias_tau_CF_k_0.05 <- mean(result_k_0.05$tau_CF-TE_k_0.05)
bias_tau_CF_k_0.1 <- mean(result_k_0.1$tau_CF-TE_k_0.1)
bias_tau_CF_k_0.2 <- mean(result_k_0.2$tau_CF-TE_k_0.2)
bias_tau_CF_k_0.3 <- mean(result_k_0.3$tau_CF-TE_k_0.3)
bias_tau_CF_k_0.4 <- mean(result_k_0.4$tau_CF-TE_k_0.4)
bias_tau_CF_k_0.7 <- mean(result_k_0.7$tau_CF-TE_k_0.7)

###G-computation estimators with CF and LASSO
bias_tau_CF_LASSO_k_0.05 <- mean(result_k_0.05$tau_CF_LASSO-TE_k_0.05)
bias_tau_CF_LASSO_k_0.1 <- mean(result_k_0.1$tau_CF_LASSO-TE_k_0.1)
bias_tau_CF_LASSO_k_0.2 <- mean(result_k_0.2$tau_CF_LASSO-TE_k_0.2)
bias_tau_CF_LASSO_k_0.3 <- mean(result_k_0.3$tau_CF_LASSO-TE_k_0.3)
bias_tau_CF_LASSO_k_0.4 <- mean(result_k_0.4$tau_CF_LASSO-TE_k_0.4)
bias_tau_CF_LASSO_k_0.7 <- mean(result_k_0.7$tau_CF_LASSO-TE_k_0.7)


###HOIF_1-based covariate adjusted estimators
bias_tau_HOIF_1_k_0.05 <- mean(result_k_0.05$tau_HOIF_1-TE_k_0.05)
bias_tau_HOIF_1_k_0.1 <- mean(result_k_0.1$tau_HOIF_1-TE_k_0.1)
bias_tau_HOIF_1_k_0.2 <- mean(result_k_0.2$tau_HOIF_1-TE_k_0.2)
bias_tau_HOIF_1_k_0.3 <- mean(result_k_0.3$tau_HOIF_1-TE_k_0.3)
bias_tau_HOIF_1_k_0.4 <- mean(result_k_0.4$tau_HOIF_1-TE_k_0.4)
bias_tau_HOIF_1_k_0.7 <- mean(result_k_0.7$tau_HOIF_1-TE_k_0.7)

###HOIF_2-based covariate adjusted estimators
bias_tau_HOIF_2_k_0.05 <- mean(result_k_0.05$tau_HOIF_2-TE_k_0.05)
bias_tau_HOIF_2_k_0.1 <- mean(result_k_0.1$tau_HOIF_2-TE_k_0.1)
bias_tau_HOIF_2_k_0.2 <- mean(result_k_0.2$tau_HOIF_2-TE_k_0.2)
bias_tau_HOIF_2_k_0.3 <- mean(result_k_0.3$tau_HOIF_2-TE_k_0.3)
bias_tau_HOIF_2_k_0.4 <- mean(result_k_0.4$tau_HOIF_2-TE_k_0.4)
bias_tau_HOIF_2_k_0.7 <- mean(result_k_0.7$tau_HOIF_2-TE_k_0.7)


### JASA estimator
bias_tau_JASA_k_0.05 <- mean(result_k_0.05$tau_JASA-TE_k_0.05)
bias_tau_JASA_k_0.1 <- mean(result_k_0.1$tau_JASA-TE_k_0.1)
bias_tau_JASA_k_0.2 <- mean(result_k_0.2$tau_JASA-TE_k_0.2)
bias_tau_JASA_k_0.3 <- mean(result_k_0.3$tau_JASA-TE_k_0.3)
bias_tau_JASA_k_0.4 <- mean(result_k_0.4$tau_JASA-TE_k_0.4)
bias_tau_JASA_k_0.7 <- mean(result_k_0.7$tau_JASA-TE_k_0.7)

### JASA_cal estimator
bias_tau_JASA_cal_k_0.05 <- mean(result_k_0.05$tau_JASA_cal-TE_k_0.05)
bias_tau_JASA_cal_k_0.1 <- mean(result_k_0.1$tau_JASA_cal-TE_k_0.1)
bias_tau_JASA_cal_k_0.2 <- mean(result_k_0.2$tau_JASA_cal-TE_k_0.2)
bias_tau_JASA_cal_k_0.3 <- mean(result_k_0.3$tau_JASA_cal-TE_k_0.3)
bias_tau_JASA_cal_k_0.4 <- mean(result_k_0.4$tau_JASA_cal-TE_k_0.4)
bias_tau_JASA_cal_k_0.7 <- mean(result_k_0.7$tau_JASA_cal-TE_k_0.7)

### Firth estimator
bias_tau_Firth_k_0.05 <- mean(result_k_0.05$tau_Firth-TE_k_0.05)
bias_tau_Firth_k_0.1 <- mean(result_k_0.1$tau_Firth-TE_k_0.1)
bias_tau_Firth_k_0.2 <- mean(result_k_0.2$tau_Firth-TE_k_0.2)
bias_tau_Firth_k_0.3 <- mean(result_k_0.3$tau_Firth-TE_k_0.3)
bias_tau_Firth_k_0.4 <- mean(result_k_0.4$tau_Firth-TE_k_0.4)
bias_tau_Firth_k_0.7 <- mean(result_k_0.7$tau_Firth-TE_k_0.7)

### Firth_target estimator
bias_tau_Firth_target_k_0.05 <- mean(result_k_0.05$tau_Firth_target-TE_k_0.05)
bias_tau_Firth_target_k_0.1 <- mean(result_k_0.1$tau_Firth_target-TE_k_0.1)
bias_tau_Firth_target_k_0.2 <- mean(result_k_0.2$tau_Firth_target-TE_k_0.2)
bias_tau_Firth_target_k_0.3 <- mean(result_k_0.3$tau_Firth_target-TE_k_0.3)
bias_tau_Firth_target_k_0.4 <- mean(result_k_0.4$tau_Firth_target-TE_k_0.4)
bias_tau_Firth_target_k_0.7 <- mean(result_k_0.7$tau_Firth_target-TE_k_0.7)


###MSE
#Unadjusted estimators
MSE_tau_unadj_k_0.05 <- mean((result_k_0.05$tau_unadj-TE_k_0.05)^2)
MSE_tau_unadj_k_0.1 <- mean((result_k_0.1$tau_unadj-TE_k_0.1)^2)
MSE_tau_unadj_k_0.2 <- mean((result_k_0.2$tau_unadj-TE_k_0.2)^2)
MSE_tau_unadj_k_0.3 <- mean((result_k_0.3$tau_unadj-TE_k_0.3)^2)
MSE_tau_unadj_k_0.4 <- mean((result_k_0.4$tau_unadj-TE_k_0.4)^2)
MSE_tau_unadj_k_0.7 <- mean((result_k_0.7$tau_unadj-TE_k_0.7)^2)


##Standard G-computation estimators
MSE_tau_stand_k_0.05 <- mean((result_k_0.05$tau_stand-TE_k_0.05)^2)
MSE_tau_stand_k_0.1 <- mean((result_k_0.1$tau_stand-TE_k_0.1)^2)
MSE_tau_stand_k_0.2 <- mean((result_k_0.2$tau_stand-TE_k_0.2)^2)
MSE_tau_stand_k_0.3 <- mean((result_k_0.3$tau_stand-TE_k_0.3)^2)
MSE_tau_stand_k_0.4 <- mean((result_k_0.4$tau_stand-TE_k_0.4)^2)
MSE_tau_stand_k_0.7 <- mean((result_k_0.7$tau_stand-TE_k_0.7)^2)

###G-computation estimators with LASSO
MSE_tau_LASSO_k_0.05 <- mean((result_k_0.05$tau_LASSO-TE_k_0.05)^2)
MSE_tau_LASSO_k_0.1 <- mean((result_k_0.1$tau_LASSO-TE_k_0.1)^2)
MSE_tau_LASSO_k_0.2 <- mean((result_k_0.2$tau_LASSO-TE_k_0.2)^2)
MSE_tau_LASSO_k_0.3 <- mean((result_k_0.3$tau_LASSO-TE_k_0.3)^2)
MSE_tau_LASSO_k_0.4 <- mean((result_k_0.4$tau_LASSO-TE_k_0.4)^2)
MSE_tau_LASSO_k_0.7 <- mean((result_k_0.7$tau_LASSO-TE_k_0.7)^2)

###G-computation estimators with CF
MSE_tau_CF_k_0.05 <- mean((result_k_0.05$tau_CF-TE_k_0.05)^2)
MSE_tau_CF_k_0.1 <- mean((result_k_0.1$tau_CF-TE_k_0.1)^2)
MSE_tau_CF_k_0.2 <- mean((result_k_0.2$tau_CF-TE_k_0.2)^2)
MSE_tau_CF_k_0.3 <- mean((result_k_0.3$tau_CF-TE_k_0.3)^2)
MSE_tau_CF_k_0.4 <- mean((result_k_0.4$tau_CF-TE_k_0.4)^2)
MSE_tau_CF_k_0.7 <- mean((result_k_0.7$tau_CF-TE_k_0.7)^2)

###G-computation estimators with CF and LASSO
MSE_tau_CF_LASSO_k_0.05 <- mean((result_k_0.05$tau_CF_LASSO-TE_k_0.05)^2)
MSE_tau_CF_LASSO_k_0.1 <- mean((result_k_0.1$tau_CF_LASSO-TE_k_0.1)^2)
MSE_tau_CF_LASSO_k_0.2 <- mean((result_k_0.2$tau_CF_LASSO-TE_k_0.2)^2)
MSE_tau_CF_LASSO_k_0.3 <- mean((result_k_0.3$tau_CF_LASSO-TE_k_0.3)^2)
MSE_tau_CF_LASSO_k_0.4 <- mean((result_k_0.4$tau_CF_LASSO-TE_k_0.4)^2)
MSE_tau_CF_LASSO_k_0.7 <- mean((result_k_0.7$tau_CF_LASSO-TE_k_0.7)^2)


###HOIF_1-based covariate adjusted estimators
MSE_tau_HOIF_1_k_0.05 <- mean((result_k_0.05$tau_HOIF_1-TE_k_0.05)^2)
MSE_tau_HOIF_1_k_0.1 <- mean((result_k_0.1$tau_HOIF_1-TE_k_0.1)^2)
MSE_tau_HOIF_1_k_0.2 <- mean((result_k_0.2$tau_HOIF_1-TE_k_0.2)^2)
MSE_tau_HOIF_1_k_0.3 <- mean((result_k_0.3$tau_HOIF_1-TE_k_0.3)^2)
MSE_tau_HOIF_1_k_0.4 <- mean((result_k_0.4$tau_HOIF_1-TE_k_0.4)^2)
MSE_tau_HOIF_1_k_0.7 <- mean((result_k_0.7$tau_HOIF_1-TE_k_0.7)^2)


###HOIF_2-based covariate adjusted estimators
MSE_tau_HOIF_2_k_0.05 <- mean((result_k_0.05$tau_HOIF_2-TE_k_0.05)^2)
MSE_tau_HOIF_2_k_0.1 <- mean((result_k_0.1$tau_HOIF_2-TE_k_0.1)^2)
MSE_tau_HOIF_2_k_0.2 <- mean((result_k_0.2$tau_HOIF_2-TE_k_0.2)^2)
MSE_tau_HOIF_2_k_0.3 <- mean((result_k_0.3$tau_HOIF_2-TE_k_0.3)^2)
MSE_tau_HOIF_2_k_0.4 <- mean((result_k_0.4$tau_HOIF_2-TE_k_0.4)^2)
MSE_tau_HOIF_2_k_0.7 <- mean((result_k_0.7$tau_HOIF_2-TE_k_0.7)^2)


### JASA estimator
MSE_tau_JASA_k_0.05 <- mean((result_k_0.05$tau_JASA-TE_k_0.05)^2)
MSE_tau_JASA_k_0.1 <- mean((result_k_0.1$tau_JASA-TE_k_0.1)^2)
MSE_tau_JASA_k_0.2 <- mean((result_k_0.2$tau_JASA-TE_k_0.2)^2)
MSE_tau_JASA_k_0.3 <- mean((result_k_0.3$tau_JASA-TE_k_0.3)^2)
MSE_tau_JASA_k_0.4 <- mean((result_k_0.4$tau_JASA-TE_k_0.4)^2)
MSE_tau_JASA_k_0.7 <- mean((result_k_0.7$tau_JASA-TE_k_0.7)^2)

### JASA_cal estimator
MSE_tau_JASA_cal_k_0.05 <- mean((result_k_0.05$tau_JASA_cal-TE_k_0.05)^2)
MSE_tau_JASA_cal_k_0.1 <- mean((result_k_0.1$tau_JASA_cal-TE_k_0.1)^2)
MSE_tau_JASA_cal_k_0.2 <- mean((result_k_0.2$tau_JASA_cal-TE_k_0.2)^2)
MSE_tau_JASA_cal_k_0.3 <- mean((result_k_0.3$tau_JASA_cal-TE_k_0.3)^2)
MSE_tau_JASA_cal_k_0.4 <- mean((result_k_0.4$tau_JASA_cal-TE_k_0.4)^2)
MSE_tau_JASA_cal_k_0.7 <- mean((result_k_0.7$tau_JASA_cal-TE_k_0.7)^2)

### Firth estimator
MSE_tau_Firth_k_0.05 <- mean((result_k_0.05$tau_Firth-TE_k_0.05)^2)
MSE_tau_Firth_k_0.1 <- mean((result_k_0.1$tau_Firth-TE_k_0.1)^2)
MSE_tau_Firth_k_0.2 <- mean((result_k_0.2$tau_Firth-TE_k_0.2)^2)
MSE_tau_Firth_k_0.3 <- mean((result_k_0.3$tau_Firth-TE_k_0.3)^2)
MSE_tau_Firth_k_0.4 <- mean((result_k_0.4$tau_Firth-TE_k_0.4)^2)
MSE_tau_Firth_k_0.7 <- mean((result_k_0.7$tau_Firth-TE_k_0.7)^2)

### Firth_target estimator
MSE_tau_Firth_target_k_0.05 <- mean((result_k_0.05$tau_Firth_target-TE_k_0.05)^2)
MSE_tau_Firth_target_k_0.1 <- mean((result_k_0.1$tau_Firth_target-TE_k_0.1)^2)
MSE_tau_Firth_target_k_0.2 <- mean((result_k_0.2$tau_Firth_target-TE_k_0.2)^2)
MSE_tau_Firth_target_k_0.3 <- mean((result_k_0.3$tau_Firth_target-TE_k_0.3)^2)
MSE_tau_Firth_target_k_0.4 <- mean((result_k_0.4$tau_Firth_target-TE_k_0.4)^2)
MSE_tau_Firth_target_k_0.7 <- mean((result_k_0.7$tau_Firth_target-TE_k_0.7)^2)


##Monte Carlo standard deviation
##Unadjusted estimators
SD_tau_unadj_k_0.05 <- sqrt(var(result_k_0.05$tau_unadj))
SD_tau_unadj_k_0.1 <- sqrt(var(result_k_0.1$tau_unadj))
SD_tau_unadj_k_0.2 <- sqrt(var(result_k_0.2$tau_unadj))
SD_tau_unadj_k_0.3 <- sqrt(var(result_k_0.3$tau_unadj))
SD_tau_unadj_k_0.4 <- sqrt(var(result_k_0.4$tau_unadj))
SD_tau_unadj_k_0.7 <- sqrt(var(result_k_0.7$tau_unadj))


###Standard G-computation estimators
SD_tau_stand_k_0.05 <- sqrt(var(result_k_0.05$tau_stand))
SD_tau_stand_k_0.1 <- sqrt(var(result_k_0.1$tau_stand))
SD_tau_stand_k_0.2 <- sqrt(var(result_k_0.2$tau_stand))
SD_tau_stand_k_0.3 <- sqrt(var(result_k_0.3$tau_stand))
SD_tau_stand_k_0.4 <- sqrt(var(result_k_0.4$tau_stand))
SD_tau_stand_k_0.7 <- sqrt(var(result_k_0.7$tau_stand))

###G-computation estimators with LASSO
SD_tau_LASSO_k_0.05 <- sqrt(var(result_k_0.05$tau_LASSO))
SD_tau_LASSO_k_0.1 <- sqrt(var(result_k_0.1$tau_LASSO))
SD_tau_LASSO_k_0.2 <- sqrt(var(result_k_0.2$tau_LASSO))
SD_tau_LASSO_k_0.3 <- sqrt(var(result_k_0.3$tau_LASSO))
SD_tau_LASSO_k_0.4 <- sqrt(var(result_k_0.4$tau_LASSO))
SD_tau_LASSO_k_0.7 <- sqrt(var(result_k_0.7$tau_LASSO))

###G-computation estimators with CF
SD_tau_CF_k_0.05 <- sqrt(var(result_k_0.05$tau_CF))
SD_tau_CF_k_0.1 <- sqrt(var(result_k_0.1$tau_CF))
SD_tau_CF_k_0.2 <- sqrt(var(result_k_0.2$tau_CF))
SD_tau_CF_k_0.3 <- sqrt(var(result_k_0.3$tau_CF))
SD_tau_CF_k_0.4 <- sqrt(var(result_k_0.4$tau_CF))
SD_tau_CF_k_0.7 <- sqrt(var(result_k_0.7$tau_CF))

###G-computation estimators with CF and LASSO
SD_tau_CF_LASSO_k_0.05 <- sqrt(var(result_k_0.05$tau_CF_LASSO))
SD_tau_CF_LASSO_k_0.1 <- sqrt(var(result_k_0.1$tau_CF_LASSO))
SD_tau_CF_LASSO_k_0.2 <- sqrt(var(result_k_0.2$tau_CF_LASSO))
SD_tau_CF_LASSO_k_0.3 <- sqrt(var(result_k_0.3$tau_CF_LASSO))
SD_tau_CF_LASSO_k_0.4 <- sqrt(var(result_k_0.4$tau_CF_LASSO))
SD_tau_CF_LASSO_k_0.7 <- sqrt(var(result_k_0.7$tau_CF_LASSO))


###HOIF_1-based covariate adjusted estimators
SD_tau_HOIF_1_k_0.05 <- sqrt(var(result_k_0.05$tau_HOIF_1))
SD_tau_HOIF_1_k_0.1 <- sqrt(var(result_k_0.1$tau_HOIF_1))
SD_tau_HOIF_1_k_0.2 <- sqrt(var(result_k_0.2$tau_HOIF_1))
SD_tau_HOIF_1_k_0.3 <- sqrt(var(result_k_0.3$tau_HOIF_1))
SD_tau_HOIF_1_k_0.4 <- sqrt(var(result_k_0.4$tau_HOIF_1))
SD_tau_HOIF_1_k_0.7 <- sqrt(var(result_k_0.7$tau_HOIF_1))


###HOIF_2-based covariate adjusted estimators
SD_tau_HOIF_2_k_0.05 <- sqrt(var(result_k_0.05$tau_HOIF_2))
SD_tau_HOIF_2_k_0.1 <- sqrt(var(result_k_0.1$tau_HOIF_2))
SD_tau_HOIF_2_k_0.2 <- sqrt(var(result_k_0.2$tau_HOIF_2))
SD_tau_HOIF_2_k_0.3 <- sqrt(var(result_k_0.3$tau_HOIF_2))
SD_tau_HOIF_2_k_0.4 <- sqrt(var(result_k_0.4$tau_HOIF_2))
SD_tau_HOIF_2_k_0.7 <- sqrt(var(result_k_0.7$tau_HOIF_2))


### JASA estimator
SD_tau_JASA_k_0.05 <- sqrt(var(result_k_0.05$tau_JASA))
SD_tau_JASA_k_0.1 <- sqrt(var(result_k_0.1$tau_JASA))
SD_tau_JASA_k_0.2 <- sqrt(var(result_k_0.2$tau_JASA))
SD_tau_JASA_k_0.3 <- sqrt(var(result_k_0.3$tau_JASA))
SD_tau_JASA_k_0.4 <- sqrt(var(result_k_0.4$tau_JASA))
SD_tau_JASA_k_0.7 <- sqrt(var(result_k_0.7$tau_JASA))

### JASA_cal estimator
SD_tau_JASA_cal_k_0.05 <- sqrt(var(result_k_0.05$tau_JASA_cal))
SD_tau_JASA_cal_k_0.1 <- sqrt(var(result_k_0.1$tau_JASA_cal))
SD_tau_JASA_cal_k_0.2 <- sqrt(var(result_k_0.2$tau_JASA_cal))
SD_tau_JASA_cal_k_0.3 <- sqrt(var(result_k_0.3$tau_JASA_cal))
SD_tau_JASA_cal_k_0.4 <- sqrt(var(result_k_0.4$tau_JASA_cal))
SD_tau_JASA_cal_k_0.7 <- sqrt(var(result_k_0.7$tau_JASA_cal))

### Firth estimator
SD_tau_Firth_k_0.05 <- sqrt(var(result_k_0.05$tau_Firth))
SD_tau_Firth_k_0.1 <- sqrt(var(result_k_0.1$tau_Firth))
SD_tau_Firth_k_0.2 <- sqrt(var(result_k_0.2$tau_Firth))
SD_tau_Firth_k_0.3 <- sqrt(var(result_k_0.3$tau_Firth))
SD_tau_Firth_k_0.4 <- sqrt(var(result_k_0.4$tau_Firth))
SD_tau_Firth_k_0.7 <- sqrt(var(result_k_0.7$tau_Firth))

### Firth_target estimator
SD_tau_Firth_target_k_0.05 <- sqrt(var(result_k_0.05$tau_Firth_target))
SD_tau_Firth_target_k_0.1 <- sqrt(var(result_k_0.1$tau_Firth_target))
SD_tau_Firth_target_k_0.2 <- sqrt(var(result_k_0.2$tau_Firth_target))
SD_tau_Firth_target_k_0.3 <- sqrt(var(result_k_0.3$tau_Firth_target))
SD_tau_Firth_target_k_0.4 <- sqrt(var(result_k_0.4$tau_Firth_target))
SD_tau_Firth_target_k_0.7 <- sqrt(var(result_k_0.7$tau_Firth_target))


###Estimated standard errors
#Unadjusted
#Based on the standard influence function
se_ave_tau_unadj_k_0.05 <- mean(result_k_0.05$se_tau_unadj_inf)
se_ave_tau_unadj_k_0.1 <- mean(result_k_0.1$se_tau_unadj_inf)
se_ave_tau_unadj_k_0.2 <- mean(result_k_0.2$se_tau_unadj_inf)
se_ave_tau_unadj_k_0.3 <- mean(result_k_0.3$se_tau_unadj_inf)
se_ave_tau_unadj_k_0.4 <- mean(result_k_0.4$se_tau_unadj_inf)
se_ave_tau_unadj_k_0.7 <- mean(result_k_0.7$se_tau_unadj_inf)

#Based on small sample correction
se_ave_tau_unadj_k_0.05_inf_corr <- mean(result_k_0.05$se_tau_unadj_inf_corr)
se_ave_tau_unadj_k_0.1_inf_corr <- mean(result_k_0.1$se_tau_unadj_inf_corr)
se_ave_tau_unadj_k_0.2_inf_corr <- mean(result_k_0.2$se_tau_unadj_inf_corr)
se_ave_tau_unadj_k_0.3_inf_corr <- mean(result_k_0.3$se_tau_unadj_inf_corr)
se_ave_tau_unadj_k_0.4_inf_corr <- mean(result_k_0.4$se_tau_unadj_inf_corr)
se_ave_tau_unadj_k_0.7_inf_corr <- mean(result_k_0.7$se_tau_unadj_inf_corr)



#Standard G-computation estimators
#Based on the standard influence function
se_ave_tau_stand_k_0.05 <- mean(result_k_0.05$se_tau_stand_inf)
se_ave_tau_stand_k_0.1 <- mean(result_k_0.1$se_tau_stand_inf)
se_ave_tau_stand_k_0.2 <- mean(result_k_0.2$se_tau_stand_inf)
se_ave_tau_stand_k_0.3 <- mean(result_k_0.3$se_tau_stand_inf)
se_ave_tau_stand_k_0.4 <- mean(result_k_0.4$se_tau_stand_inf)
se_ave_tau_stand_k_0.7 <- mean(result_k_0.7$se_tau_stand_inf)

#Based on small sample correction
se_ave_tau_stand_k_0.05_inf_corr <- mean(result_k_0.05$se_tau_stand_inf_corr)
se_ave_tau_stand_k_0.1_inf_corr <- mean(result_k_0.1$se_tau_stand_inf_corr)
se_ave_tau_stand_k_0.2_inf_corr <- mean(result_k_0.2$se_tau_stand_inf_corr)
se_ave_tau_stand_k_0.3_inf_corr <- mean(result_k_0.3$se_tau_stand_inf_corr)
se_ave_tau_stand_k_0.4_inf_corr <- mean(result_k_0.4$se_tau_stand_inf_corr)
se_ave_tau_stand_k_0.7_inf_corr <- mean(result_k_0.7$se_tau_stand_inf_corr)


###G-computation estimators with LASSO
#Based on the standard influence function
se_ave_tau_LASSO_k_0.05 <- mean(result_k_0.05$se_tau_LASSO_inf)
se_ave_tau_LASSO_k_0.1 <- mean(result_k_0.1$se_tau_LASSO_inf)
se_ave_tau_LASSO_k_0.2 <- mean(result_k_0.2$se_tau_LASSO_inf)
se_ave_tau_LASSO_k_0.3 <- mean(result_k_0.3$se_tau_LASSO_inf)
se_ave_tau_LASSO_k_0.4 <- mean(result_k_0.4$se_tau_LASSO_inf)
se_ave_tau_LASSO_k_0.7 <- mean(result_k_0.7$se_tau_LASSO_inf)

#Based on small sample correction
se_ave_tau_LASSO_k_0.05_inf_corr <- mean(result_k_0.05$se_tau_LASSO_inf_corr)
se_ave_tau_LASSO_k_0.1_inf_corr <- mean(result_k_0.1$se_tau_LASSO_inf_corr)
se_ave_tau_LASSO_k_0.2_inf_corr <- mean(result_k_0.2$se_tau_LASSO_inf_corr)
se_ave_tau_LASSO_k_0.3_inf_corr <- mean(result_k_0.3$se_tau_LASSO_inf_corr)
se_ave_tau_LASSO_k_0.4_inf_corr <- mean(result_k_0.4$se_tau_LASSO_inf_corr)
se_ave_tau_LASSO_k_0.7_inf_corr <- mean(result_k_0.7$se_tau_LASSO_inf_corr)


###G-computation estimators with CF
#Based on the standard influence function
se_ave_tau_CF_k_0.05 <- mean(result_k_0.05$se_tau_CF_inf)
se_ave_tau_CF_k_0.1 <- mean(result_k_0.1$se_tau_CF_inf)
se_ave_tau_CF_k_0.2 <- mean(result_k_0.2$se_tau_CF_inf)
se_ave_tau_CF_k_0.3 <- mean(result_k_0.3$se_tau_CF_inf)
se_ave_tau_CF_k_0.4 <- mean(result_k_0.4$se_tau_CF_inf)
se_ave_tau_CF_k_0.7 <- mean(result_k_0.7$se_tau_CF_inf)

#Based on small sample correction
se_ave_tau_CF_k_0.05_inf_corr <- mean(result_k_0.05$se_tau_CF_inf_corr)
se_ave_tau_CF_k_0.1_inf_corr <- mean(result_k_0.1$se_tau_CF_inf_corr)
se_ave_tau_CF_k_0.2_inf_corr <- mean(result_k_0.2$se_tau_CF_inf_corr)
se_ave_tau_CF_k_0.3_inf_corr <- mean(result_k_0.3$se_tau_CF_inf_corr)
se_ave_tau_CF_k_0.4_inf_corr <- mean(result_k_0.4$se_tau_CF_inf_corr)
se_ave_tau_CF_k_0.7_inf_corr <- mean(result_k_0.7$se_tau_CF_inf_corr)

###G-computation estimators with CF and LASSO
#Based on the standard influence function
se_ave_tau_CF_LASSO_k_0.05 <- mean(result_k_0.05$se_tau_CF_LASSO_inf)
se_ave_tau_CF_LASSO_k_0.1 <- mean(result_k_0.1$se_tau_CF_LASSO_inf)
se_ave_tau_CF_LASSO_k_0.2 <- mean(result_k_0.2$se_tau_CF_LASSO_inf)
se_ave_tau_CF_LASSO_k_0.3 <- mean(result_k_0.3$se_tau_CF_LASSO_inf)
se_ave_tau_CF_LASSO_k_0.4 <- mean(result_k_0.4$se_tau_CF_LASSO_inf)
se_ave_tau_CF_LASSO_k_0.7 <- mean(result_k_0.7$se_tau_CF_LASSO_inf)


###HOIF_1-based covariate adjusted estimators
#Based on the standard influence function
se_ave_tau_HOIF_1_k_0.05 <- mean(result_k_0.05$se_tau_HOIF_1)
se_ave_tau_HOIF_1_k_0.1 <- mean(result_k_0.1$se_tau_HOIF_1)
se_ave_tau_HOIF_1_k_0.2 <- mean(result_k_0.2$se_tau_HOIF_1)
se_ave_tau_HOIF_1_k_0.3 <- mean(result_k_0.3$se_tau_HOIF_1)
se_ave_tau_HOIF_1_k_0.4 <- mean(result_k_0.4$se_tau_HOIF_1)
se_ave_tau_HOIF_1_k_0.7 <- mean(result_k_0.7$se_tau_HOIF_1)


###HOIF_2-based covariate adjusted estimators
#Based on the standard influence function
se_ave_tau_HOIF_2_k_0.05 <- mean(result_k_0.05$se_tau_HOIF_2)
se_ave_tau_HOIF_2_k_0.1 <- mean(result_k_0.1$se_tau_HOIF_2)
se_ave_tau_HOIF_2_k_0.2 <- mean(result_k_0.2$se_tau_HOIF_2)
se_ave_tau_HOIF_2_k_0.3 <- mean(result_k_0.3$se_tau_HOIF_2)
se_ave_tau_HOIF_2_k_0.4 <- mean(result_k_0.4$se_tau_HOIF_2)
se_ave_tau_HOIF_2_k_0.7 <- mean(result_k_0.7$se_tau_HOIF_2)


### JASA estimator
se_ave_tau_JASA_k_0.05 <- mean(result_k_0.05$se_tau_JASA)
se_ave_tau_JASA_k_0.1 <- mean(result_k_0.1$se_tau_JASA)
se_ave_tau_JASA_k_0.2 <- mean(result_k_0.2$se_tau_JASA)
se_ave_tau_JASA_k_0.3 <- mean(result_k_0.3$se_tau_JASA)
se_ave_tau_JASA_k_0.4 <- mean(result_k_0.4$se_tau_JASA)
se_ave_tau_JASA_k_0.7 <- mean(result_k_0.7$se_tau_JASA)

### JASA_cal estimator
se_ave_tau_JASA_cal_k_0.05 <- mean(result_k_0.05$se_tau_JASA_cal)
se_ave_tau_JASA_cal_k_0.1 <- mean(result_k_0.1$se_tau_JASA_cal)
se_ave_tau_JASA_cal_k_0.2 <- mean(result_k_0.2$se_tau_JASA_cal)
se_ave_tau_JASA_cal_k_0.3 <- mean(result_k_0.3$se_tau_JASA_cal)
se_ave_tau_JASA_cal_k_0.4 <- mean(result_k_0.4$se_tau_JASA_cal)
se_ave_tau_JASA_cal_k_0.7 <- mean(result_k_0.7$se_tau_JASA_cal)

### Firth estimator
se_ave_tau_Firth_k_0.05 <- mean(result_k_0.05$se_tau_Firth)
se_ave_tau_Firth_k_0.1 <- mean(result_k_0.1$se_tau_Firth)
se_ave_tau_Firth_k_0.2 <- mean(result_k_0.2$se_tau_Firth)
se_ave_tau_Firth_k_0.3 <- mean(result_k_0.3$se_tau_Firth)
se_ave_tau_Firth_k_0.4 <- mean(result_k_0.4$se_tau_Firth)
se_ave_tau_Firth_k_0.7 <- mean(result_k_0.7$se_tau_Firth)

### Firth_target estimator
se_ave_tau_Firth_target_k_0.05 <- mean(result_k_0.05$se_tau_Firth_target)
se_ave_tau_Firth_target_k_0.1 <- mean(result_k_0.1$se_tau_Firth_target)
se_ave_tau_Firth_target_k_0.2 <- mean(result_k_0.2$se_tau_Firth_target)
se_ave_tau_Firth_target_k_0.3 <- mean(result_k_0.3$se_tau_Firth_target)
se_ave_tau_Firth_target_k_0.4 <- mean(result_k_0.4$se_tau_Firth_target)
se_ave_tau_Firth_target_k_0.7 <- mean(result_k_0.7$se_tau_Firth_target)


##Average width of confidence intervals
#The variance estimator is based on the influence function
##Unadjusted estimators
ME <- 2*1.96
WCI_tau_unadj_k_0.05 <- mean(ME*result_k_0.05$se_tau_unadj_inf)
WCI_tau_unadj_k_0.1 <- mean(ME*result_k_0.1$se_tau_unadj_inf)
WCI_tau_unadj_k_0.2 <- mean(ME*result_k_0.2$se_tau_unadj_inf)
WCI_tau_unadj_k_0.3 <- mean(ME*result_k_0.3$se_tau_unadj_inf)
WCI_tau_unadj_k_0.4 <- mean(ME*result_k_0.4$se_tau_unadj_inf)
WCI_tau_unadj_k_0.7 <- mean(ME*result_k_0.7$se_tau_unadj_inf)


###Standard G-computation estimators
WCI_tau_stand_k_0.05 <- mean(ME*result_k_0.05$se_tau_stand_inf)
WCI_tau_stand_k_0.1 <- mean(ME*result_k_0.1$se_tau_stand_inf)
WCI_tau_stand_k_0.2 <- mean(ME*result_k_0.2$se_tau_stand_inf)
WCI_tau_stand_k_0.3 <- mean(ME*result_k_0.3$se_tau_stand_inf)
WCI_tau_stand_k_0.4 <- mean(ME*result_k_0.4$se_tau_stand_inf)
WCI_tau_stand_k_0.7 <- mean(ME*result_k_0.7$se_tau_stand_inf)


###G-computation estimators with LASSO
WCI_tau_LASSO_k_0.05 <- mean(ME*result_k_0.05$se_tau_LASSO_inf)
WCI_tau_LASSO_k_0.1 <- mean(ME*result_k_0.1$se_tau_LASSO_inf)
WCI_tau_LASSO_k_0.2 <- mean(ME*result_k_0.2$se_tau_LASSO_inf)
WCI_tau_LASSO_k_0.3 <- mean(ME*result_k_0.3$se_tau_LASSO_inf)
WCI_tau_LASSO_k_0.4 <- mean(ME*result_k_0.4$se_tau_LASSO_inf)
WCI_tau_LASSO_k_0.7 <- mean(ME*result_k_0.7$se_tau_LASSO_inf)


###G-computation estimators with CF
WCI_tau_CF_k_0.05 <- mean(ME*result_k_0.05$se_tau_CF_inf)
WCI_tau_CF_k_0.1 <- mean(ME*result_k_0.1$se_tau_CF_inf)
WCI_tau_CF_k_0.2 <- mean(ME*result_k_0.2$se_tau_CF_inf)
WCI_tau_CF_k_0.3 <- mean(ME*result_k_0.3$se_tau_CF_inf)
WCI_tau_CF_k_0.4 <- mean(ME*result_k_0.4$se_tau_CF_inf)
WCI_tau_CF_k_0.7 <- mean(ME*result_k_0.7$se_tau_CF_inf)


###G-computation estimators with CF and LASSO
WCI_tau_CF_LASSO_k_0.05 <- mean(ME*result_k_0.05$se_tau_CF_LASSO_inf)
WCI_tau_CF_LASSO_k_0.1 <- mean(ME*result_k_0.1$se_tau_CF_LASSO_inf)
WCI_tau_CF_LASSO_k_0.2 <- mean(ME*result_k_0.2$se_tau_CF_LASSO_inf)
WCI_tau_CF_LASSO_k_0.3 <- mean(ME*result_k_0.3$se_tau_CF_LASSO_inf)
WCI_tau_CF_LASSO_k_0.4 <- mean(ME*result_k_0.4$se_tau_CF_LASSO_inf)
WCI_tau_CF_LASSO_k_0.7 <- mean(ME*result_k_0.7$se_tau_CF_LASSO_inf)


###HOIF_1-motivated covariate adjusted estimators
WCI_tau_HOIF_1_k_0.05 <- mean(ME*result_k_0.05$se_tau_HOIF_1)
WCI_tau_HOIF_1_k_0.1 <- mean(ME*result_k_0.1$se_tau_HOIF_1)
WCI_tau_HOIF_1_k_0.2 <- mean(ME*result_k_0.2$se_tau_HOIF_1)
WCI_tau_HOIF_1_k_0.3 <- mean(ME*result_k_0.3$se_tau_HOIF_1)
WCI_tau_HOIF_1_k_0.4 <- mean(ME*result_k_0.4$se_tau_HOIF_1)
WCI_tau_HOIF_1_k_0.7 <- mean(ME*result_k_0.7$se_tau_HOIF_1)


###HOIF_2-motivated covariate adjusted estimators
WCI_tau_HOIF_2_k_0.05 <- mean(ME*result_k_0.05$se_tau_HOIF_2)
WCI_tau_HOIF_2_k_0.1 <- mean(ME*result_k_0.1$se_tau_HOIF_2)
WCI_tau_HOIF_2_k_0.2 <- mean(ME*result_k_0.2$se_tau_HOIF_2)
WCI_tau_HOIF_2_k_0.3 <- mean(ME*result_k_0.3$se_tau_HOIF_2)
WCI_tau_HOIF_2_k_0.4 <- mean(ME*result_k_0.4$se_tau_HOIF_2)
WCI_tau_HOIF_2_k_0.7 <- mean(ME*result_k_0.7$se_tau_HOIF_2)

### JASA estimator
WCI_tau_JASA_k_0.05 <- mean(ME*result_k_0.05$se_tau_JASA)
WCI_tau_JASA_k_0.1 <- mean(ME*result_k_0.1$se_tau_JASA)
WCI_tau_JASA_k_0.2 <- mean(ME*result_k_0.2$se_tau_JASA)
WCI_tau_JASA_k_0.3 <- mean(ME*result_k_0.3$se_tau_JASA)
WCI_tau_JASA_k_0.4 <- mean(ME*result_k_0.4$se_tau_JASA)
WCI_tau_JASA_k_0.7 <- mean(ME*result_k_0.7$se_tau_JASA)

### JASA_cal estimator
WCI_tau_JASA_cal_k_0.05 <- mean(ME*result_k_0.05$se_tau_JASA_cal)
WCI_tau_JASA_cal_k_0.1 <- mean(ME*result_k_0.1$se_tau_JASA_cal)
WCI_tau_JASA_cal_k_0.2 <- mean(ME*result_k_0.2$se_tau_JASA_cal)
WCI_tau_JASA_cal_k_0.3 <- mean(ME*result_k_0.3$se_tau_JASA_cal)
WCI_tau_JASA_cal_k_0.4 <- mean(ME*result_k_0.4$se_tau_JASA_cal)
WCI_tau_JASA_cal_k_0.7 <- mean(ME*result_k_0.7$se_tau_JASA_cal)

### Firth estimator
WCI_tau_Firth_k_0.05 <- mean(ME*result_k_0.05$se_tau_Firth)
WCI_tau_Firth_k_0.1 <- mean(ME*result_k_0.1$se_tau_Firth)
WCI_tau_Firth_k_0.2 <- mean(ME*result_k_0.2$se_tau_Firth)
WCI_tau_Firth_k_0.3 <- mean(ME*result_k_0.3$se_tau_Firth)
WCI_tau_Firth_k_0.4 <- mean(ME*result_k_0.4$se_tau_Firth)
WCI_tau_Firth_k_0.7 <- mean(ME*result_k_0.7$se_tau_Firth)

### Firth_target estimator
WCI_tau_Firth_target_k_0.05 <- mean(ME*result_k_0.05$se_tau_Firth_target)
WCI_tau_Firth_target_k_0.1 <- mean(ME*result_k_0.1$se_tau_Firth_target)
WCI_tau_Firth_target_k_0.2 <- mean(ME*result_k_0.2$se_tau_Firth_target)
WCI_tau_Firth_target_k_0.3 <- mean(ME*result_k_0.3$se_tau_Firth_target)
WCI_tau_Firth_target_k_0.4 <- mean(ME*result_k_0.4$se_tau_Firth_target)
WCI_tau_Firth_target_k_0.7 <- mean(ME*result_k_0.7$se_tau_Firth_target)




##Average width of confidence intervals
#The variance estimator is based on the influence function with small-sample correction
###Unadjusted G-computation estimators
WCI_tau_unadj_k_0.05_inf_corr <- mean(ME*result_k_0.05$se_tau_unadj_inf_corr)
WCI_tau_unadj_k_0.1_inf_corr <- mean(ME*result_k_0.1$se_tau_unadj_inf_corr)
WCI_tau_unadj_k_0.2_inf_corr <- mean(ME*result_k_0.2$se_tau_unadj_inf_corr)
WCI_tau_unadj_k_0.3_inf_corr <- mean(ME*result_k_0.3$se_tau_unadj_inf_corr)
WCI_tau_unadj_k_0.4_inf_corr <- mean(ME*result_k_0.4$se_tau_unadj_inf_corr)
WCI_tau_unadj_k_0.7_inf_corr <- mean(ME*result_k_0.7$se_tau_unadj_inf_corr)

###Standard G-computation estimators
WCI_tau_stand_k_0.05_inf_corr <- mean(ME*result_k_0.05$se_tau_stand_inf_corr)
WCI_tau_stand_k_0.1_inf_corr <- mean(ME*result_k_0.1$se_tau_stand_inf_corr)
WCI_tau_stand_k_0.2_inf_corr <- mean(ME*result_k_0.2$se_tau_stand_inf_corr)
WCI_tau_stand_k_0.3_inf_corr <- mean(ME*result_k_0.3$se_tau_stand_inf_corr)
WCI_tau_stand_k_0.4_inf_corr <- mean(ME*result_k_0.4$se_tau_stand_inf_corr)
WCI_tau_stand_k_0.7_inf_corr <- mean(ME*result_k_0.7$se_tau_stand_inf_corr)


###G-computation estimators with LASSO
WCI_tau_LASSO_k_0.05_inf_corr <- mean(ME*result_k_0.05$se_tau_LASSO_inf_corr)
WCI_tau_LASSO_k_0.1_inf_corr <- mean(ME*result_k_0.1$se_tau_LASSO_inf_corr)
WCI_tau_LASSO_k_0.2_inf_corr <- mean(ME*result_k_0.2$se_tau_LASSO_inf_corr)
WCI_tau_LASSO_k_0.3_inf_corr <- mean(ME*result_k_0.3$se_tau_LASSO_inf_corr)
WCI_tau_LASSO_k_0.4_inf_corr <- mean(ME*result_k_0.4$se_tau_LASSO_inf_corr)
WCI_tau_LASSO_k_0.7_inf_corr <- mean(ME*result_k_0.7$se_tau_LASSO_inf_corr)


###G-computation estimators with CF
WCI_tau_CF_k_0.05_inf_corr <- mean(ME*result_k_0.05$se_tau_CF_inf_corr)
WCI_tau_CF_k_0.1_inf_corr <- mean(ME*result_k_0.1$se_tau_CF_inf_corr)
WCI_tau_CF_k_0.2_inf_corr <- mean(ME*result_k_0.2$se_tau_CF_inf_corr)
WCI_tau_CF_k_0.3_inf_corr <- mean(ME*result_k_0.3$se_tau_CF_inf_corr)
WCI_tau_CF_k_0.4_inf_corr <- mean(ME*result_k_0.4$se_tau_CF_inf_corr)
WCI_tau_CF_k_0.7_inf_corr <- mean(ME*result_k_0.7$se_tau_CF_inf_corr)


###Coverage function
#Function to calculate coverage probability
true_parameter_in_ci<-function(data,point_estimate,se,true,nrep){
  true_in_ci<-sapply(1:nrep,function(i){
    lower_ci<-data[i,point_estimate]-1.96*data[i,se]
    upper_ci<-data[i,point_estimate]+1.96*data[i,se]
    if(true>=lower_ci & true<=upper_ci){1}else{0}
  })
  cov<-100*sum(true_in_ci)/nrep
  return(list("cov"=cov))
}


#Unadjusted estimators
#Based on the standard influence function
cov_tau_unadj_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_unadj_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_unadj_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_unadj_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_unadj_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_unadj_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_unadj","se_tau_unadj_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

###Standard G-computation estimators
#Based on the standard influence function
cov_tau_stand_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_stand","se_tau_stand_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_stand_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_stand","se_tau_stand_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_stand_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_stand","se_tau_stand_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_stand_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_stand","se_tau_stand_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_stand_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_stand","se_tau_stand_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_stand_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_stand","se_tau_stand_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

#Based on small-sample correction
cov_tau_stand_k_0.05_inf_corr<-true_parameter_in_ci(result_k_0.05,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_stand_k_0.1_inf_corr<-true_parameter_in_ci(result_k_0.1,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_stand_k_0.2_inf_corr<-true_parameter_in_ci(result_k_0.2,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_stand_k_0.3_inf_corr<-true_parameter_in_ci(result_k_0.3,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_stand_k_0.4_inf_corr<-true_parameter_in_ci(result_k_0.4,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_stand_k_0.7_inf_corr<-true_parameter_in_ci(result_k_0.7,"tau_stand","se_tau_stand_inf_corr",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]


###G-computation estimators with LASSO
#Based on the standard influence function
cov_tau_LASSO_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_LASSO_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_LASSO_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_LASSO_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_LASSO_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_LASSO_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_LASSO","se_tau_LASSO_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

#Based on small-sample correction
cov_tau_LASSO_k_0.05_inf_corr<-true_parameter_in_ci(result_k_0.05,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_LASSO_k_0.1_inf_corr<-true_parameter_in_ci(result_k_0.1,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_LASSO_k_0.2_inf_corr<-true_parameter_in_ci(result_k_0.2,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_LASSO_k_0.3_inf_corr<-true_parameter_in_ci(result_k_0.3,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_LASSO_k_0.4_inf_corr<-true_parameter_in_ci(result_k_0.4,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_LASSO_k_0.7_inf_corr<-true_parameter_in_ci(result_k_0.7,"tau_LASSO","se_tau_LASSO_inf_corr",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]


###G-computation estimators with CF
#Based on the standard influence function
cov_tau_CF_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_CF","se_tau_CF_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_CF_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_CF","se_tau_CF_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_CF_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_CF","se_tau_CF_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_CF_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_CF","se_tau_CF_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_CF_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_CF","se_tau_CF_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_CF_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_CF","se_tau_CF_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

#Based on small-sample correction
cov_tau_CF_k_0.05_inf_corr<-true_parameter_in_ci(result_k_0.05,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_CF_k_0.1_inf_corr<-true_parameter_in_ci(result_k_0.1,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_CF_k_0.2_inf_corr<-true_parameter_in_ci(result_k_0.2,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_CF_k_0.3_inf_corr<-true_parameter_in_ci(result_k_0.3,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_CF_k_0.4_inf_corr<-true_parameter_in_ci(result_k_0.4,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_CF_k_0.7_inf_corr<-true_parameter_in_ci(result_k_0.7,"tau_CF","se_tau_CF_inf_corr",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

###G-computation estimators with CF and LASSO
#Based on the standard influence function
cov_tau_CF_LASSO_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_CF_LASSO_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_CF_LASSO_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_CF_LASSO_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_CF_LASSO_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_CF_LASSO_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]



###G-computation estimators with CF and LASSO
#Based on the standard influence function
cov_tau_CF_LASSO_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_CF_LASSO_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_CF_LASSO_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_CF_LASSO_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_CF_LASSO_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_CF_LASSO_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_CF_LASSO","se_tau_CF_LASSO_inf",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

###HOIF_1-based covariate adjusted estimators
#Based on the standard influence function
cov_tau_HOIF_1_k_0.05<-true_parameter_in_ci(result_k_0.05,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_HOIF_1_k_0.1<-true_parameter_in_ci(result_k_0.1,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_HOIF_1_k_0.2<-true_parameter_in_ci(result_k_0.2,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_HOIF_1_k_0.3<-true_parameter_in_ci(result_k_0.3,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_HOIF_1_k_0.4<-true_parameter_in_ci(result_k_0.4,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_HOIF_1_k_0.7<-true_parameter_in_ci(result_k_0.7,"tau_HOIF_1","se_tau_HOIF_1",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

###HOIF_2-based covariate adjusted estimators
cov_tau_HOIF_2_k_0.05<-true_parameter_in_ci(result_k_0.05,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_HOIF_2_k_0.1<-true_parameter_in_ci(result_k_0.1,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_HOIF_2_k_0.2<-true_parameter_in_ci(result_k_0.2,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_HOIF_2_k_0.3<-true_parameter_in_ci(result_k_0.3,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_HOIF_2_k_0.4<-true_parameter_in_ci(result_k_0.4,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_HOIF_2_k_0.7<-true_parameter_in_ci(result_k_0.7,"tau_HOIF_2","se_tau_HOIF_2",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]
#######################################################################################

### JASA estimator
cov_tau_JASA_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_JASA","se_tau_JASA",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_JASA_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_JASA","se_tau_JASA",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_JASA_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_JASA","se_tau_JASA",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_JASA_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_JASA","se_tau_JASA",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_JASA_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_JASA","se_tau_JASA",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_JASA_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_JASA","se_tau_JASA",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

### JASA_cal estimator
cov_tau_JASA_cal_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_JASA_cal_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_JASA_cal_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_JASA_cal_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_JASA_cal_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_JASA_cal_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_JASA_cal","se_tau_JASA_cal",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

### Firth estimator
cov_tau_Firth_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_Firth","se_tau_Firth",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_Firth_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_Firth","se_tau_Firth",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_Firth_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_Firth","se_tau_Firth",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_Firth_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_Firth","se_tau_Firth",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_Firth_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_Firth","se_tau_Firth",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_Firth_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_Firth","se_tau_Firth",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]

### Firth_target estimator
cov_tau_Firth_target_k_0.05 <-true_parameter_in_ci(result_k_0.05,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.05,dim(result_k_0.05)[1])[[1]]
cov_tau_Firth_target_k_0.1 <-true_parameter_in_ci(result_k_0.1,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.1,dim(result_k_0.1)[1])[[1]]
cov_tau_Firth_target_k_0.2 <-true_parameter_in_ci(result_k_0.2,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.2,dim(result_k_0.2)[1])[[1]]
cov_tau_Firth_target_k_0.3 <-true_parameter_in_ci(result_k_0.3,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.3,dim(result_k_0.3)[1])[[1]]
cov_tau_Firth_target_k_0.4 <-true_parameter_in_ci(result_k_0.4,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.4,dim(result_k_0.4)[1])[[1]]
cov_tau_Firth_target_k_0.7 <-true_parameter_in_ci(result_k_0.7,"tau_Firth_target","se_tau_Firth_target",true=TE_k_0.7,dim(result_k_0.7)[1])[[1]]


###Tables
#Bias and MSE
Y1_k_0.05=0.6820663;Y1_k_0.1=0.6756631;Y1_k_0.2=0.6666502;Y1_k_0.3=0.6621297;Y1_k_0.4=0.6594322;Y1_k_0.7=0.65341
Y0_k_0.05=0.5470243;Y0_k_0.1=0.5457495;Y0_k_0.2=0.5429026;Y0_k_0.3=0.5417;Y0_k_0.4=0.5413291;Y0_k_0.7=0.5393423

df_bias_MSE <- data.frame(  
  k=c("0.05","0.1","0.2","0.3","0.4","0.7"),
  Y1=c(Y1_k_0.05,Y1_k_0.1,Y1_k_0.2,Y1_k_0.3,Y1_k_0.4,Y1_k_0.7),
  Y0=c(Y0_k_0.05,Y0_k_0.1,Y0_k_0.2,Y0_k_0.3,Y0_k_0.4,Y0_k_0.7),
  bias_tau_unadj=c(bias_tau_unadj_k_0.05,bias_tau_unadj_k_0.1,bias_tau_unadj_k_0.2,bias_tau_unadj_k_0.3,bias_tau_unadj_k_0.4,bias_tau_unadj_k_0.7),
  bias_tau_stand=c(bias_tau_stand_k_0.05,bias_tau_stand_k_0.1,bias_tau_stand_k_0.2,bias_tau_stand_k_0.3,bias_tau_stand_k_0.4,bias_tau_stand_k_0.7),
  bias_tau_LASSO=c(bias_tau_LASSO_k_0.05,bias_tau_LASSO_k_0.1,bias_tau_LASSO_k_0.2,bias_tau_LASSO_k_0.3,bias_tau_LASSO_k_0.4,bias_tau_LASSO_k_0.7),
  bias_tau_CF=c(bias_tau_CF_k_0.05,bias_tau_CF_k_0.1,bias_tau_CF_k_0.2,bias_tau_CF_k_0.3,bias_tau_CF_k_0.4,bias_tau_CF_k_0.7),
  bias_tau_CF_LASSO=c(bias_tau_CF_LASSO_k_0.05,bias_tau_CF_LASSO_k_0.1,bias_tau_CF_LASSO_k_0.2,bias_tau_CF_LASSO_k_0.3,bias_tau_CF_LASSO_k_0.4,bias_tau_CF_LASSO_k_0.7),
  bias_tau_HOIF_1=c(bias_tau_HOIF_1_k_0.05,bias_tau_HOIF_1_k_0.1,bias_tau_HOIF_1_k_0.2,bias_tau_HOIF_1_k_0.3,bias_tau_HOIF_1_k_0.4,bias_tau_HOIF_1_k_0.7),
  bias_tau_HOIF_2=c(bias_tau_HOIF_2_k_0.05,bias_tau_HOIF_2_k_0.1,bias_tau_HOIF_2_k_0.2,bias_tau_HOIF_2_k_0.3,bias_tau_HOIF_2_k_0.4,bias_tau_HOIF_2_k_0.7),
  bias_tau_JASA=c(bias_tau_JASA_k_0.05,bias_tau_JASA_k_0.1,bias_tau_JASA_k_0.2,bias_tau_JASA_k_0.3,bias_tau_JASA_k_0.4,bias_tau_JASA_k_0.7),
  bias_tau_JASA_cal=c(bias_tau_JASA_cal_k_0.05,bias_tau_JASA_cal_k_0.1,bias_tau_JASA_cal_k_0.2,bias_tau_JASA_cal_k_0.3,bias_tau_JASA_cal_k_0.4,bias_tau_JASA_cal_k_0.7),
  bias_tau_Firth=c(bias_tau_Firth_k_0.05,bias_tau_Firth_k_0.1,bias_tau_Firth_k_0.2,bias_tau_Firth_k_0.3,bias_tau_Firth_k_0.4,bias_tau_Firth_k_0.7),
  bias_tau_Firth_target=c(bias_tau_Firth_target_k_0.05,bias_tau_Firth_target_k_0.1,bias_tau_Firth_target_k_0.2,bias_tau_Firth_target_k_0.3,bias_tau_Firth_target_k_0.4,bias_tau_Firth_target_k_0.7),
  
  
  
  
  MSE_tau_unadj=c(MSE_tau_unadj_k_0.05,MSE_tau_unadj_k_0.1,MSE_tau_unadj_k_0.2,MSE_tau_unadj_k_0.3,MSE_tau_unadj_k_0.4,MSE_tau_unadj_k_0.7),
  MSE_tau_stand=c(MSE_tau_stand_k_0.05,MSE_tau_stand_k_0.1,MSE_tau_stand_k_0.2,MSE_tau_stand_k_0.3,MSE_tau_stand_k_0.4,MSE_tau_stand_k_0.7),
  MSE_tau_LASSO=c(MSE_tau_LASSO_k_0.05,MSE_tau_LASSO_k_0.1,MSE_tau_LASSO_k_0.2,MSE_tau_LASSO_k_0.3,MSE_tau_LASSO_k_0.4,MSE_tau_LASSO_k_0.7),
  MSE_tau_CF=c(MSE_tau_CF_k_0.05,MSE_tau_CF_k_0.1,MSE_tau_CF_k_0.2,MSE_tau_CF_k_0.3,MSE_tau_CF_k_0.4,MSE_tau_CF_k_0.7),
  MSE_tau_CF_LASSO=c(MSE_tau_CF_LASSO_k_0.05,MSE_tau_CF_LASSO_k_0.1,MSE_tau_CF_LASSO_k_0.2,MSE_tau_CF_LASSO_k_0.3,MSE_tau_CF_LASSO_k_0.4,MSE_tau_CF_LASSO_k_0.7),
  MSE_tau_HOIF_1=c(MSE_tau_HOIF_1_k_0.05,MSE_tau_HOIF_1_k_0.1,MSE_tau_HOIF_1_k_0.2,MSE_tau_HOIF_1_k_0.3,MSE_tau_HOIF_1_k_0.4,MSE_tau_HOIF_1_k_0.7),
  MSE_tau_HOIF_2=c(MSE_tau_HOIF_2_k_0.05,MSE_tau_HOIF_2_k_0.1,MSE_tau_HOIF_2_k_0.2,MSE_tau_HOIF_2_k_0.3,MSE_tau_HOIF_2_k_0.4,MSE_tau_HOIF_2_k_0.7),
  MSE_tau_JASA=c(MSE_tau_JASA_k_0.05,MSE_tau_JASA_k_0.1,MSE_tau_JASA_k_0.2,MSE_tau_JASA_k_0.3,MSE_tau_JASA_k_0.4,MSE_tau_JASA_k_0.7),
  MSE_tau_JASA_cal=c(MSE_tau_JASA_cal_k_0.05,MSE_tau_JASA_cal_k_0.1,MSE_tau_JASA_cal_k_0.2,MSE_tau_JASA_cal_k_0.3,MSE_tau_JASA_cal_k_0.4,MSE_tau_JASA_cal_k_0.7),
  MSE_tau_Firth=c(MSE_tau_Firth_k_0.05,MSE_tau_Firth_k_0.1,MSE_tau_Firth_k_0.2,MSE_tau_Firth_k_0.3,MSE_tau_Firth_k_0.4,MSE_tau_Firth_k_0.7),
  MSE_tau_Firth_target=c(MSE_tau_Firth_target_k_0.05,MSE_tau_Firth_target_k_0.1,MSE_tau_Firth_target_k_0.2,MSE_tau_Firth_target_k_0.3,MSE_tau_Firth_target_k_0.4,MSE_tau_Firth_target_k_0.7))

#Convert the data frame to LaTeX format
latex_table_bias_MSE <- xtable(df_bias_MSE,digits=3)

# Print the LaTeX table
print(latex_table_bias_MSE,include.rownames = FALSE)
##########################################################################################################

#Monte Carlo standard deviation and estimated standard error based on influence function
df_SD_se_inf <- data.frame(  
  k=c("0.05","0.1","0.2","0.3","0.4","0.7"),
  SD_tau_unadj=c(SD_tau_unadj_k_0.05,SD_tau_unadj_k_0.1,SD_tau_unadj_k_0.2,SD_tau_unadj_k_0.3,SD_tau_unadj_k_0.4,SD_tau_unadj_k_0.7),
  SD_tau_stand=c(SD_tau_stand_k_0.05,SD_tau_stand_k_0.1,SD_tau_stand_k_0.2,SD_tau_stand_k_0.3,SD_tau_stand_k_0.4,SD_tau_stand_k_0.7),
  SD_tau_LASSO=c(SD_tau_LASSO_k_0.05,SD_tau_LASSO_k_0.1,SD_tau_LASSO_k_0.2,SD_tau_LASSO_k_0.3,SD_tau_LASSO_k_0.4,SD_tau_LASSO_k_0.7),
  SD_tau_CF=c(SD_tau_CF_k_0.05,SD_tau_CF_k_0.1,SD_tau_CF_k_0.2,SD_tau_CF_k_0.3,SD_tau_CF_k_0.4,SD_tau_CF_k_0.7),
  SD_tau_CF_LASSO=c(SD_tau_CF_LASSO_k_0.05,SD_tau_CF_LASSO_k_0.1,SD_tau_CF_LASSO_k_0.2,SD_tau_CF_LASSO_k_0.3,SD_tau_CF_LASSO_k_0.4,SD_tau_CF_LASSO_k_0.7),
  SD_tau_HOIF_1=c(SD_tau_HOIF_1_k_0.05,SD_tau_HOIF_1_k_0.1,SD_tau_HOIF_1_k_0.2,SD_tau_HOIF_1_k_0.3,SD_tau_HOIF_1_k_0.4,SD_tau_HOIF_1_k_0.7),
  SD_tau_HOIF_2=c(SD_tau_HOIF_2_k_0.05,SD_tau_HOIF_2_k_0.1,SD_tau_HOIF_2_k_0.2,SD_tau_HOIF_2_k_0.3,SD_tau_HOIF_2_k_0.4,SD_tau_HOIF_2_k_0.7),
  SD_tau_JASA=c(SD_tau_JASA_k_0.05,SD_tau_JASA_k_0.1,SD_tau_JASA_k_0.2,SD_tau_JASA_k_0.3,SD_tau_JASA_k_0.4,SD_tau_JASA_k_0.7),
  SD_tau_JASA_cal=c(SD_tau_JASA_cal_k_0.05,SD_tau_JASA_cal_k_0.1,SD_tau_JASA_cal_k_0.2,SD_tau_JASA_cal_k_0.3,SD_tau_JASA_cal_k_0.4,SD_tau_JASA_cal_k_0.7),
  SD_tau_Firth=c(SD_tau_Firth_k_0.05,SD_tau_Firth_k_0.1,SD_tau_Firth_k_0.2,SD_tau_Firth_k_0.3,SD_tau_Firth_k_0.4,SD_tau_Firth_k_0.7),
  SD_tau_Firth_target=c(SD_tau_Firth_target_k_0.05,SD_tau_Firth_target_k_0.1,SD_tau_Firth_target_k_0.2,SD_tau_Firth_target_k_0.3,SD_tau_Firth_target_k_0.4,SD_tau_Firth_target_k_0.7),
  
  
  
  se_tau_unadj=c(se_ave_tau_unadj_k_0.05,se_ave_tau_unadj_k_0.1,se_ave_tau_unadj_k_0.2,se_ave_tau_unadj_k_0.3,se_ave_tau_unadj_k_0.4,se_ave_tau_unadj_k_0.7),
  se_inf_tau_stand=c(se_ave_tau_stand_k_0.05,se_ave_tau_stand_k_0.1,se_ave_tau_stand_k_0.2,se_ave_tau_stand_k_0.3,se_ave_tau_stand_k_0.4,se_ave_tau_stand_k_0.7),
  se_inf_tau_LASSO=c(se_ave_tau_LASSO_k_0.05,se_ave_tau_LASSO_k_0.1,se_ave_tau_LASSO_k_0.2,se_ave_tau_LASSO_k_0.3,se_ave_tau_LASSO_k_0.4,se_ave_tau_LASSO_k_0.7),
  se_inf_tau_CF=c(se_ave_tau_CF_k_0.05,se_ave_tau_CF_k_0.1,se_ave_tau_CF_k_0.2,se_ave_tau_CF_k_0.3,se_ave_tau_CF_k_0.4,se_ave_tau_CF_k_0.7),
  se_inf_tau_CF_LASSO=c(se_ave_tau_CF_LASSO_k_0.05,se_ave_tau_CF_LASSO_k_0.1,se_ave_tau_CF_LASSO_k_0.2,se_ave_tau_CF_LASSO_k_0.3,se_ave_tau_CF_LASSO_k_0.4,se_ave_tau_CF_LASSO_k_0.7),
  se_inf_tau_HOIF_1=c(se_ave_tau_HOIF_1_k_0.05,se_ave_tau_HOIF_1_k_0.1,se_ave_tau_HOIF_1_k_0.2,se_ave_tau_HOIF_1_k_0.3,se_ave_tau_HOIF_1_k_0.4,se_ave_tau_HOIF_1_k_0.7),
  se_inf_tau_HOIF_2=c(se_ave_tau_HOIF_2_k_0.05,se_ave_tau_HOIF_2_k_0.1,se_ave_tau_HOIF_2_k_0.2,se_ave_tau_HOIF_2_k_0.3,se_ave_tau_HOIF_2_k_0.4,se_ave_tau_HOIF_2_k_0.7),
  se_tau_JASA=c(se_ave_tau_JASA_k_0.05,se_ave_tau_JASA_k_0.1,se_ave_tau_JASA_k_0.2,se_ave_tau_JASA_k_0.3,se_ave_tau_JASA_k_0.4,se_ave_tau_JASA_k_0.7),
  se_tau_JASA_cal=c(se_ave_tau_JASA_cal_k_0.05,se_ave_tau_JASA_cal_k_0.1,se_ave_tau_JASA_cal_k_0.2,se_ave_tau_JASA_cal_k_0.3,se_ave_tau_JASA_cal_k_0.4,se_ave_tau_JASA_cal_k_0.7),
  se_tau_Firth=c(se_ave_tau_Firth_k_0.05,se_ave_tau_Firth_k_0.1,se_ave_tau_Firth_k_0.2,se_ave_tau_Firth_k_0.3,se_ave_tau_Firth_k_0.4,se_ave_tau_Firth_k_0.7),
  se_tau_Firth_target=c(se_ave_tau_Firth_target_k_0.05,se_ave_tau_Firth_target_k_0.1,se_ave_tau_Firth_target_k_0.2,se_ave_tau_Firth_target_k_0.3,se_ave_tau_Firth_target_k_0.4,se_ave_tau_Firth_target_k_0.7))

#Convert the data frame to LaTeX format
latex_table_SD_se_inf <- xtable(df_SD_se_inf,digits=3)

# Print the LaTeX table
print(latex_table_SD_se_inf,include.rownames = FALSE)
##########################################################################################################
#Average width of CIs based the standard error estimated using the influence function and Coverage probabilities
df_WCI_inf_cov <- data.frame(  
  k=c("0.05","0.1","0.2","0.3","0.4","0.7"),
  WCI_tau_unadj=c(WCI_tau_unadj_k_0.05,WCI_tau_unadj_k_0.1,WCI_tau_unadj_k_0.2,WCI_tau_unadj_k_0.3,WCI_tau_unadj_k_0.4,WCI_tau_unadj_k_0.7),
  WCI_inf_tau_stand=c(WCI_tau_stand_k_0.05,WCI_tau_stand_k_0.1,WCI_tau_stand_k_0.2,WCI_tau_stand_k_0.3,WCI_tau_stand_k_0.4,WCI_tau_stand_k_0.7),
  WCI_inf_tau_LASSO=c(WCI_tau_LASSO_k_0.05,WCI_tau_LASSO_k_0.1,WCI_tau_LASSO_k_0.2,WCI_tau_LASSO_k_0.3,WCI_tau_LASSO_k_0.4,WCI_tau_LASSO_k_0.7),
  WCI_inf_tau_CF=c(WCI_tau_CF_k_0.05,WCI_tau_CF_k_0.1,WCI_tau_CF_k_0.2,WCI_tau_CF_k_0.3,WCI_tau_CF_k_0.4,WCI_tau_CF_k_0.7),
  WCI_inf_tau_CF_LASSO=c(WCI_tau_CF_LASSO_k_0.05,WCI_tau_CF_LASSO_k_0.1,WCI_tau_CF_LASSO_k_0.2,WCI_tau_CF_LASSO_k_0.3,WCI_tau_CF_LASSO_k_0.4,WCI_tau_CF_LASSO_k_0.7),
  WCI_inf_tau_HOIF_1=c(WCI_tau_HOIF_1_k_0.05,WCI_tau_HOIF_1_k_0.1,WCI_tau_HOIF_1_k_0.2,WCI_tau_HOIF_1_k_0.3,WCI_tau_HOIF_1_k_0.4,WCI_tau_HOIF_1_k_0.7),
  WCI_inf_tau_HOIF_2=c(WCI_tau_HOIF_2_k_0.05,WCI_tau_HOIF_2_k_0.1,WCI_tau_HOIF_2_k_0.2,WCI_tau_HOIF_2_k_0.3,WCI_tau_HOIF_2_k_0.4,WCI_tau_HOIF_2_k_0.7),
  WCI_tau_JASA=c(WCI_tau_JASA_k_0.05,WCI_tau_JASA_k_0.1,WCI_tau_JASA_k_0.2,WCI_tau_JASA_k_0.3,WCI_tau_JASA_k_0.4,WCI_tau_JASA_k_0.7),
  WCI_tau_JASA_cal=c(WCI_tau_JASA_cal_k_0.05,WCI_tau_JASA_cal_k_0.1,WCI_tau_JASA_cal_k_0.2,WCI_tau_JASA_cal_k_0.3,WCI_tau_JASA_cal_k_0.4,WCI_tau_JASA_cal_k_0.7),
  WCI_tau_Firth=c(WCI_tau_Firth_k_0.05,WCI_tau_Firth_k_0.1,WCI_tau_Firth_k_0.2,WCI_tau_Firth_k_0.3,WCI_tau_Firth_k_0.4,WCI_tau_Firth_k_0.7),
  WCI_tau_Firth_target=c(WCI_tau_Firth_target_k_0.05,WCI_tau_Firth_target_k_0.1,WCI_tau_Firth_target_k_0.2,WCI_tau_Firth_target_k_0.3,WCI_tau_Firth_target_k_0.4,WCI_tau_Firth_target_k_0.7),
  
  
  cov_inf_tau_unadj=c(cov_tau_unadj_k_0.05,cov_tau_unadj_k_0.1,cov_tau_unadj_k_0.2,cov_tau_unadj_k_0.3,cov_tau_unadj_k_0.4,cov_tau_unadj_k_0.7),
  cov_inf_tau_stand=c(cov_tau_stand_k_0.05,cov_tau_stand_k_0.1,cov_tau_stand_k_0.2,cov_tau_stand_k_0.3,cov_tau_stand_k_0.4,cov_tau_stand_k_0.7),
  cov_inf_tau_LASSO=c(cov_tau_LASSO_k_0.05,cov_tau_LASSO_k_0.1,cov_tau_LASSO_k_0.2,cov_tau_LASSO_k_0.3,cov_tau_LASSO_k_0.4,cov_tau_LASSO_k_0.7),
  cov_inf_tau_CF=c(cov_tau_CF_k_0.05,cov_tau_CF_k_0.1,cov_tau_CF_k_0.2,cov_tau_CF_k_0.3,cov_tau_CF_k_0.4,cov_tau_CF_k_0.7),
  cov_inf_tau_CF_LASSO=c(cov_tau_CF_LASSO_k_0.05,cov_tau_CF_LASSO_k_0.1,cov_tau_CF_LASSO_k_0.2,cov_tau_CF_LASSO_k_0.3,cov_tau_CF_LASSO_k_0.4,cov_tau_CF_LASSO_k_0.7),
  cov_inf_tau_HOIF_1=c(cov_tau_HOIF_1_k_0.05,cov_tau_HOIF_1_k_0.1,cov_tau_HOIF_1_k_0.2,cov_tau_HOIF_1_k_0.3,cov_tau_HOIF_1_k_0.4,cov_tau_HOIF_1_k_0.7),
  cov_inf_tau_HOIF_2=c(cov_tau_HOIF_2_k_0.05,cov_tau_HOIF_2_k_0.1,cov_tau_HOIF_2_k_0.2,cov_tau_HOIF_2_k_0.3,cov_tau_HOIF_2_k_0.4,cov_tau_HOIF_2_k_0.7),
  cov_tau_JASA=c(cov_tau_JASA_k_0.05,cov_tau_JASA_k_0.1,cov_tau_JASA_k_0.2,cov_tau_JASA_k_0.3,cov_tau_JASA_k_0.4,cov_tau_JASA_k_0.7),
  cov_tau_JASA_cal=c(cov_tau_JASA_cal_k_0.05,cov_tau_JASA_cal_k_0.1,cov_tau_JASA_cal_k_0.2,cov_tau_JASA_cal_k_0.3,cov_tau_JASA_cal_k_0.4,cov_tau_JASA_cal_k_0.7),
  cov_tau_Firth=c(cov_tau_Firth_k_0.05,cov_tau_Firth_k_0.1,cov_tau_Firth_k_0.2,cov_tau_Firth_k_0.3,cov_tau_Firth_k_0.4,cov_tau_Firth_k_0.7),
  cov_tau_Firth_target=c(cov_tau_Firth_target_k_0.05,cov_tau_Firth_target_k_0.1,cov_tau_Firth_target_k_0.2,cov_tau_Firth_target_k_0.3,cov_tau_Firth_target_k_0.4,cov_tau_Firth_target_k_0.7))


#Convert the data frame to LaTeX format
latex_table_WCI_inf_cov <- xtable(df_WCI_inf_cov,digits=3)

# Print the LaTeX table
print(latex_table_WCI_inf_cov,include.rownames = FALSE)
###########################################################################################################
#Monte Carlo standard deviation and estimated standard error based on influence function that uses degree of freedom correction
df_SD_se_inf_corr <- data.frame(  
  k=c("0.05","0.1","0.2","0.3","0.4","0.7"),
  SD_tau_stand=c(SD_tau_stand_k_0.05,SD_tau_stand_k_0.1,SD_tau_stand_k_0.2,SD_tau_stand_k_0.3,SD_tau_stand_k_0.4,SD_tau_stand_k_0.7),
  SD_tau_LASSO=c(SD_tau_LASSO_k_0.05,SD_tau_LASSO_k_0.1,SD_tau_LASSO_k_0.2,SD_tau_LASSO_k_0.3,SD_tau_LASSO_k_0.4,SD_tau_LASSO_k_0.7),
  SD_tau_CF=c(SD_tau_CF_k_0.05,SD_tau_CF_k_0.1,SD_tau_CF_k_0.2,SD_tau_CF_k_0.3,SD_tau_CF_k_0.4,SD_tau_CF_k_0.7),
  
  se_inf_corr_tau_stand=c(se_ave_tau_stand_k_0.05_inf_corr,se_ave_tau_stand_k_0.1_inf_corr,se_ave_tau_stand_k_0.2_inf_corr,se_ave_tau_stand_k_0.3_inf_corr,se_ave_tau_stand_k_0.4_inf_corr,se_ave_tau_stand_k_0.7_inf_corr),
  se_inf_corr_tau_LASSO=c(se_ave_tau_LASSO_k_0.05_inf_corr,se_ave_tau_LASSO_k_0.1_inf_corr,se_ave_tau_LASSO_k_0.2_inf_corr,se_ave_tau_LASSO_k_0.3_inf_corr,se_ave_tau_LASSO_k_0.4_inf_corr,se_ave_tau_LASSO_k_0.7_inf_corr),
  se_inf_corr_tau_CF=c(se_ave_tau_CF_k_0.05_inf_corr,se_ave_tau_CF_k_0.1_inf_corr,se_ave_tau_CF_k_0.2_inf_corr,se_ave_tau_CF_k_0.3_inf_corr,se_ave_tau_CF_k_0.4_inf_corr,se_ave_tau_CF_k_0.7_inf_corr))


#Convert the data frame to LaTeX format
latex_table_SD_se_inf_corr <- xtable(df_SD_se_inf_corr,digits=3)

# Print the LaTeX table
print(latex_table_SD_se_inf_corr,include.rownames = FALSE)
###########################################################################################################
#Average width of CIs based the standard error estimated using the influence function with df correction and Coverage probabilities
df_WCI_inf_corr_cov <- data.frame(  
  k=c("0.05","0.1","0.2","0.3","0.4","0.7"),
  WCI_inf_corr_tau_stand=c(WCI_tau_stand_k_0.05_inf_corr,WCI_tau_stand_k_0.1_inf_corr,WCI_tau_stand_k_0.2_inf_corr,WCI_tau_stand_k_0.3_inf_corr,WCI_tau_stand_k_0.4_inf_corr,WCI_tau_stand_k_0.7_inf_corr),
  WCI_inf_corr_tau_LASSO=c(WCI_tau_LASSO_k_0.05_inf_corr,WCI_tau_LASSO_k_0.1_inf_corr,WCI_tau_LASSO_k_0.2_inf_corr,WCI_tau_LASSO_k_0.3_inf_corr,WCI_tau_LASSO_k_0.4_inf_corr,WCI_tau_LASSO_k_0.7_inf_corr),
  WCI_inf_corr_tau_CF=c(WCI_tau_CF_k_0.05_inf_corr,WCI_tau_CF_k_0.1_inf_corr,WCI_tau_CF_k_0.2_inf_corr,WCI_tau_CF_k_0.3_inf_corr,WCI_tau_CF_k_0.4_inf_corr,WCI_tau_CF_k_0.7_inf_corr),
  
  cov_inf_corr_tau_stand=c(cov_tau_stand_k_0.05_inf_corr,cov_tau_stand_k_0.1_inf_corr,cov_tau_stand_k_0.2_inf_corr,cov_tau_stand_k_0.3_inf_corr,cov_tau_stand_k_0.4_inf_corr,cov_tau_stand_k_0.7_inf_corr),
  cov_inf_corr_tau_LASSO=c(cov_tau_LASSO_k_0.05_inf_corr,cov_tau_LASSO_k_0.1_inf_corr,cov_tau_LASSO_k_0.2_inf_corr,cov_tau_LASSO_k_0.3_inf_corr,cov_tau_LASSO_k_0.4_inf_corr,cov_tau_LASSO_k_0.7_inf_corr),
  cov_inf_corr_tau_CF=c(cov_tau_CF_k_0.05_inf_corr,cov_tau_CF_k_0.1_inf_corr,cov_tau_CF_k_0.2_inf_corr,cov_tau_CF_k_0.3_inf_corr,cov_tau_CF_k_0.4_inf_corr,cov_tau_CF_k_0.7_inf_corr))


#Convert the data frame to LaTeX format
latex_table_WCI_inf_corr_cov <- xtable(df_WCI_inf_corr_cov,digits=3)

# Print the LaTeX table
print(latex_table_WCI_inf_corr_cov,include.rownames = FALSE)
#********************end*****************************************************



