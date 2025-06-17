# Functions to run simulations when the nuisance generating model is linear.

even_indicator <- function(num){
  ifelse(num %% 2 == 0, 1, 0)
}

gen_nonlinear <- function(whichtype, xt, t, t_total, alpha0, lambda1){
  if (whichtype == "beta"){
    out <- alpha0 + lambda1*(dbeta(xt/6, 2, 2) + dbeta(t/t_total, 2, 2))
  } else if (whichtype == "periodic"){
    out <- alpha0 + lambda1*(sin(t) + sin(xt))
  } else if (whichtype == "steps"){
    out <- alpha0 + lambda1*(even_indicator(t) + even_indicator(floor(10*xt)))
  } else if (whichtype == "linear"){
    out <- alpha0 + lambda1*(t/t_total + xt/6)
  }
  out
}

dgm_nonlinear <- function(i_total, t_total, beta0, beta1, alpha0, lambda1, alpha1, lambda2, functype){
  dat <- c()
  for (i in 1:i_total){
    #Generate the data:
    t <- 1:t_total
    X1i <- sample(c(0, 1, 2), t_total, replace=TRUE, prob=c(0.1, 0.2, 0.7)) #time-varying covariate
    X2i <- runif(t_total, -2, 2) 
    Ai <- rbinom(t_total, 1, pt) #treatment variable
    
    delta_i <- beta0 + beta1*X2i   #the CEE
    eta_i <- gen_nonlinear(whichtype = functype, xt = X2i, t = t, t_total = t_total, 
                           alpha0 = alpha0, lambda1 = lambda1)
    miss_prob_genmod <- gen_nonlinear(whichtype = functype, xt = X2i, t = t, t_total = t_total,
                                      alpha0 = alpha1, lambda1 = lambda2)
    
    errori <- rnorm(t_total) #no correlations
    Yi <- Ai*delta_i + eta_i + errori
    
    #Missing data:
    miss_prob_true <- exp(miss_prob_genmod)/(1 + exp(miss_prob_genmod))
    Ri <- 1*rbinom(t_total, 1, miss_prob_true)
    Yi <- ifelse(Ri == 1, Yi, NA)
    
    #So the data here is {Ai, X1i, X2i, RiYi, Yi}
    dat <- rbind(dat, cbind(i, t, Ai, X1i, X2i, Ri, Yi, miss_prob_true))
  }
  colnames(dat) <- c("id", "t", "A", "X1", "X2", "R", "Y", "miss_prob_true")
  dat <- as.data.frame(dat)
  dat
}



## Estimate beta

#############################################################Method: numerical solution
beta_numerical_solution <- function(dat, moderator){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  
  ee <- function(beta){
    Sbeta <- as.numeric(S_mat %*% beta)
    
    dat$eta_hat <- (pt + pt - 1)*Sbeta + (1-pt)*dat$mu1 + pt*dat$mu0
    
    U1 <- (dat$Y - dat$eta_hat - (dat$A - pt)*Sbeta)*(dat$A - pt)
    U1 <- ifelse(is.na(U1), 0, U1)
    EU <- dat$expect_y_hat - dat$eta_hat - (dat$A - pt)*Sbeta
    
    Ustar <- rep(NA, length(beta)) # value of estimating function
    for (ibeta in 1:length(beta)) {
      Ustar[ibeta] <- sum((dat$R/dat$e_hat)*U1 * S_mat[, ibeta] 
                          - ((dat$R - dat$e_hat)/dat$e_hat)*(dat$A - pt)*S_mat[, ibeta]*EU)
    }
    
    Ustar <- Ustar/i_total
    return(Ustar)
  }
  solution <- multiroot(ee, rep(0, ncol(S_mat)),  useFortran = FALSE)
  beta_hat_solve <- solution$root
  
  beta_hat_solve
}




beta_var_calculate <- function(dat, beta, moderator, 
                               e_fit, expect_y_fit){
  
  # eta_hat for optimal formula
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  Sbeta <- as.numeric(S_mat %*% beta)
  
  dat$eta_hat <- (pt + pt - 1)*Sbeta + (1-pt)*dat$mu1 + pt*dat$mu0
  i_total <- length(unique(dat$id))
  total_person_decisionpoint <- nrow(dat)
  
  #For e, logistic regression
  X_e <- model.matrix(e_fit) # Design matrix
  mu_e <- e_fit$fitted.values # Fitted values
  W_e <- diag(mu_e * (1 - mu_e))
  #For mu, linear regression
  dat <- dat %>% 
    mutate("A:X2" = A*X2,
           "A:t" = A*t)
  variables_mu <- colnames(model.matrix(expect_y_fit))[-1]
  X_mu <- as.matrix(cbind(1, dat[, variables_mu]))
  
  p <- ncol(X_e) + ncol(X_mu) + length(beta)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  for (it in 1:total_person_decisionpoint){
    
    #calculate derivative of Utilde wrt parameter of interest and nuisance parameter
    ## w.r.t beta
    dat_it <- dat[it, ]
    fx <- matrix(c(1, dat_it[, moderator]), nrow = length(beta), ncol = 1)
    dee_Utilde_beta <- -((dat_it$A + pt - 1)*(dat_it$A - pt))*(fx%*%t(fx))  
    ## w.r.t gamma_e
    gamma_e <- e_fit$coefficients
    dee_e_gamma_e <- as.numeric(exp(-gamma_e%*%X_e[it,])/(1 + exp(-gamma_e%*%X_e[it,]))^2)*(-X_e[it,])
    dee_Utilde_gamma_e <- (dat_it$R/dat_it$e_hat^2)*(dat_it$Y - dat_it$expect_y_hat)*(dat_it$A - pt)
    dee_Utilde_gamma_e <- ifelse(is.na(dat_it$Y), 0, dee_Utilde_gamma_e)*fx%*%(dee_e_gamma_e)
    ## w.r.t gamma_mu
    gamma_mu <- expect_y_fit$coefficients
    dee_Utilde_gamma_mu <- -((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*(dat_it$A - pt)*fx%*%X_mu[it,]
    
    #calcualte derivative of nuisance parameter's own estimating equations
    ## dee of e
    dee_e <- -t(t(X_e[it,])) %*% W_e[it,it] %*% t(X_e[it,])
    ## dee of mu
    dee_mu <- -t(t(X_mu[it,])) %*% t(X_mu[it,])
    
    #combine them
    MatrixDiag <- function(A, B){
      new_rows <- nrow(A) + nrow(B)
      new_cols <- ncol(A) + ncol(B)
      # Create an empty matrix filled with zeros
      result_matrix <- matrix(0, nrow = new_rows, ncol = new_cols)
      # Place A in the top-left corner
      result_matrix[1:nrow(A), 1:ncol(A)] <- A
      # Place B in the bottom-right corner
      result_matrix[(nrow(A) + 1):new_rows, (ncol(A) + 1):new_cols] <- B
      result_matrix
    }
    ## combine them
    max_rows <- max(nrow(dee_Utilde_beta), nrow(dee_Utilde_gamma_e), nrow(dee_Utilde_gamma_mu))
    pad_matrix <- function(mat, max_rows) {
      nrow_diff <- max_rows - nrow(mat)
      if (nrow_diff > 0) {
        zero_pad <- matrix(0, nrow = nrow_diff, ncol = ncol(mat))
        mat <- rbind(mat, zero_pad)
      }
      return(mat)
    }
    dee_Utilde_beta <- pad_matrix(dee_Utilde_beta, max_rows)
    dee_Utilde_gamma_e <- pad_matrix(dee_Utilde_gamma_e, max_rows)
    dee_Utilde_gamma_mu <- pad_matrix(dee_Utilde_gamma_mu, max_rows)
    combined_matrix_bottom <- cbind(dee_Utilde_gamma_e, dee_Utilde_gamma_mu, dee_Utilde_beta)
    combined_matrix_diag <- MatrixDiag(dee_e, dee_mu)
    dee_sum[it, , ] <- rbind(cbind(combined_matrix_diag, matrix(0, 
                                                                ncol = ncol(combined_matrix_bottom)- ncol(combined_matrix_diag),
                                                                nrow = nrow(combined_matrix_diag))
    ), combined_matrix_bottom)
    
    ##########################################
    # calculate meat
    ## EE of gamma_e
    psi_gamma_e <-  (dat_it$R - dat_it$e_hat)*t(t(X_e[it,]))
    ## EE of gamma_mu
    psi_gamma_mu <- t(t(X_mu[it, ])) * ifelse(is.na(dat_it$Y), 0, dat_it$Y - dat_it$expect_y_hat)
    ## EE of beta
    u1_it <- (dat_it$Y - dat_it$eta_hat - (dat_it$A - pt)*t(fx)%*%beta)*(dat_it$A - pt)
    u1_it <- ifelse(is.na(dat_it$Y), 0, as.numeric(u1_it))*fx
    u2_it <- (dat_it$A - pt)*fx*as.numeric(dat_it$expect_y_hat - dat_it$eta_hat - (dat_it$A - pt)*t(fx)%*%beta)
    ustar_it <- (dat_it$R/dat_it$e_hat)*u1_it - ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*u2_it
    ## combine them
    meat_sum[it, ,] <- rbind(psi_gamma_e, psi_gamma_mu, ustar_it)%*%t(rbind(psi_gamma_e, psi_gamma_mu, ustar_it))
  }
  dee <- apply(dee_sum, c(2,3), sum)/i_total
  dee_inv <- solve(dee)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
  asy_var <- diag(var_cov)
  asy_var_beta <- asy_var[(p-1):p]
  
  asy_var_beta
  
}



## Estimate beta
beta_est_function <- function(dat, moderator,
                              e_formula, e_est_type, expect_y_est_type, expect_y_formula, 
                              ...){
  i_total <- length(unique(dat$id))
  ## Estimate the missingness et
  if (e_est_type == "lm"){
    e_fit <- glm(e_formula, family = binomial(link='logit'), data = dat)
    e_hat <- predict(e_fit, newdata = dat, type = "response")
    dat$e_hat <- e_hat    
  } else if (e_est_type == "gam"){
    e_fit <- gam(e_formula, method = "REML", family = binomial("logit"), data = dat)
    e_hat <- predict(e_fit, newdata = dat, type = "response")
    dat$e_hat <- e_hat   
  }
  
  
  ## Estimate E(Yt|Ht, At)
  if (expect_y_est_type == "lm"){
    expect_y_fit <- lm(expect_y_formula, data = dat)
    expect_y_hat <- predict(expect_y_fit, newdata = dat, type = "response")
    dat$expect_y_hat <- expect_y_hat
    
    newdata_for_mu1_mu0 <-  data.frame(rbind(cbind(A = 1, t = dat$t, moderator = dat[, moderator]),
                                             cbind(A = 0, t = dat$t, moderator = dat[, moderator])))
    colnames(newdata_for_mu1_mu0) <- c("A", "t", moderator)
    mu_s <- cbind(newdata_for_mu1_mu0, 
                  expect_y = predict(expect_y_fit, newdata = newdata_for_mu1_mu0, type = "response"))
    mu1_hat <- mu_s[mu_s$A == 1, "expect_y"]
    mu0_hat <- mu_s[mu_s$A == 0, "expect_y"]
    dat$mu1_hat <- mu1_hat
    dat$mu0_hat <- mu0_hat
  } else if (expect_y_est_type == "gam"){
    if (expect_y_est_gam_withA == "byA"){
      dat_a0 <- subset(dat, A == 0)
      dat_a1 <- subset(dat, A == 1)
      fit_gam_a0 <- gam(expect_y_formula, data = dat_a0)
      fit_gam_a1 <- gam(expect_y_formula, data = dat_a1)
      
      dat_set_all_a_to_0 <- dat
      dat_set_all_a_to_0$A <- 0
      mu0_hat <- predict(fit_gam_a0, newdata = dat_set_all_a_to_0, type = "response")
      
      dat_set_all_a_to_1 <- dat
      dat_set_all_a_to_1$A <- 1
      mu1_hat <- predict(fit_gam_a1, newdata = dat_set_all_a_to_1, type = "response")
      
      dat$expect_y_hat <- ifelse(dat$A, mu1_hat, mu0_hat)
      dat$mu0_hat <- mu0_hat
      dat$mu1_hat <- mu1_hat
    } else if(expect_y_est_gam_withA == "notbyA"){
      fit_gam <- gam(expect_y_formula, data = dat)
      dat_set_all_a_to_0 <- dat
      dat_set_all_a_to_0$A <- 0
      mu0_hat <- predict(fit_gam, newdata = dat_set_all_a_to_0, type = "response")
      
      dat_set_all_a_to_1 <- dat
      dat_set_all_a_to_1$A <- 1
      mu1_hat <- predict(fit_gam, newdata = dat_set_all_a_to_1, type = "response")
      
      dat$expect_y_hat <- ifelse(dat$A, mu1_hat, mu0_hat)
      dat$mu0_hat <- mu0_hat
      dat$mu1_hat <- mu1_hat
    }
  }
  
  ## Estimate beta
  beta_est_numerical <- beta_numerical_solution(dat = dat, moderator = moderator)
  
  ## Estimate the variance
  beta_var <- beta_var_calculate(dat = dat, beta = beta_est_numerical, moderator = moderator,
                                 e_fit = e_fit, expect_y_fit = expect_y_fit)
  
  beta_est <- list(beta_est_numerical, beta_var)
  names(beta_est) <- c("beta_est_numerical", "asy_var")
  beta_est
}



eval_ci <- function(beta_hat, beta_var, beta_true){
  ci_low <- beta_hat - 1.96*sqrt(beta_var)
  ci_high <- beta_hat + 1.96*sqrt(beta_var)
  beta_0_in <- ifelse(beta_true[1] <= ci_high[1] & beta_true[1] >= ci_low[1], 1, 0)
  beta_1_in <- ifelse(beta_true[2] <= ci_high[2] & beta_true[2] >= ci_low[2], 1, 0)
  beta_in <- ifelse(beta_0_in == 1 & beta_1_in == 1, 1, 0)
  list(beta_0_in, beta_1_in, beta_in)
}





# simulation function
sim_nonlineardgm <- function(i_total_list, t_total, beta0, beta1, 
                             alpha0, alpha3, alpha4, lambda1, functype,
                             moderator, 
                             expect_y_est_type, expect_y_formula, e_formula, e_est_type, ...){
  result_allsim <- list()
  result_i_total <- list()
  for (i_total in i_total_list){
    
    for (isim in 1:nsim){
      if (isim %% 10 == 0) {
        cat(isim, "")
      }
      dat <- dgm_nonlinear(i_total = i_total, t_total = t_total,  
                           beta0 = beta0, beta1 = beta1, 
                           alpha0 = alpha0, lambda1 = lambda1,
                           alpha1 = alpha1, lambda2 = lambda2,
                           functype = functype)
      
      beta_hat <- beta_est_function(dat = dat, moderator = moderator,
                                    e_formula = e_formula, 
                                    e_est_type = e_est_type,
                                    expect_y_est_type = expect_y_est_type,
                                    expect_y_formula = expect_y_formula, 
                                    expect_y_est_gam_withA = expect_y_est_gam_withA)
      
      ci_in <- eval_ci(beta_hat$beta_est_numerical, beta_hat$asy_var, beta_true)
      
      #collect the result
      result <- list(beta_hat, ci_in, mean(dat$miss_prob_true), 
                     i_total = i_total, t_total = t_total)
      result_allsim[[isim]] <- result
    }
    
    result_i_total <- c(result_i_total, list(result_allsim))
  }
  result_i_total
}

### performance evaluation
beta_result_function <- function(beta_hat_allsim, beta_true){
  beta_hat <- apply(beta_hat_allsim, 1, mean)
  bias <- apply(beta_hat_allsim - beta_true, 1, mean)
  mse <- apply((beta_hat_allsim - beta_true)^2, 1, mean)
  result <- list(beta_hat, bias, mse)
  names(result) <- c("beta_hat", "bias", "mse")
  result
}


organize_sim <- function(result_i_total, beta_true){
  result_organized <- list()
  for (i in 1:length(result_i_total)){
    result_allsim <- result_i_total[[i]]
    i_total <- result_allsim[[1]]$i_total
    t_total <- result_allsim[[1]]$t_total
    
    # beta_hat
    beta_result_allsim <-  (sapply(result_allsim, "[[", 1))
    beta_hat_allsim <- do.call(cbind, beta_result_allsim["beta_est_numerical",])
    beta_numerical_result <- beta_result_function(beta_hat_allsim, beta_true = beta_true)
    
    #asymptotic variance
    asy_var_allsim <- do.call(cbind, beta_result_allsim["asy_var",])
    se_result <- apply(sqrt(asy_var_allsim), 1, mean)
    asy_var <- apply(asy_var_allsim, 1, mean)
    beta_numerical_result$se <- se_result
    
    
    #derivative and meat
    # dee_allsim <- beta_result_allsim["dee",]
    # dee_hat <- Reduce("+", dee_allsim) / length(dee_allsim) #dee_hat
    # deeinv_allsim <- beta_result_allsim["dee_inv",]
    # deeinv_hat <- Reduce("+", deeinv_allsim) / length(deeinv_allsim) #deeinv_hat
    # meat_allsim <- beta_result_allsim["meat",]
    # meat_hat <- Reduce("+", meat_allsim) / length(meat_allsim) #meat_hat
    # var_check <- list(dee_hat = dee_hat, deeinv_hat = deeinv_hat, meat_hat = meat_hat)
    
    #CP based on estimated variance
    ci_in_allsim <- matrix(unlist(sapply(result_allsim, "[[", 2)), nrow = 3, byrow = FALSE)
    cp_result <- apply(ci_in_allsim, 1, mean)
    # CP based on true variance
    true_varaince <-  apply(beta_hat_allsim, 1, var)
    ci_truevariance <- list()
    for (i in 1:ncol(beta_hat_allsim)){
      ci_truevariance[[i]] <- eval_ci(beta_hat = beta_hat_allsim[,i],
                                      beta_var = true_varaince, 
                                      beta_true = beta_true)
    }
    cp_result_truevariance <- c(mean(sapply(ci_truevariance, "[[", 1)),
                                mean(sapply(ci_truevariance, "[[", 2)))
    
    beta_numerical_result$asy_var <- asy_var
    beta_numerical_result$true_var <- true_varaince
    
    #True observed response probability
    et_result <- mean(sapply(result_allsim, "[[", 3))
    
    #collect results
    hold <- list(i_total, t_total,
                 beta_numerical_result,
                 cp_result,
                 et_result,
                 cp_result_truevariance)
    names(hold) <- c("i_total", "t_total",
                     "beta_numerical_result",
                     "coverage_probability", 
                     "missing_prob_true",
                     "cp_result_truevariance")
    result_organized <- c(result_organized, list(hold))
  }
  result_organized
}

## Run simulation
library(rootSolve) 
library(mgcv)
library(tidyverse)
#True parameter values
beta0 <- 1.5
beta1 <- 2.1
beta_true <- c(beta0, beta1)
alpha0 <- 0.5
lambda1 <- 1.5
alpha1 <- 0.5
lambda2 <- 1.5

t_total <- 20
pt = 0.4
i_total_list <- c(50, 100, 150, 200)
nsim <- 1000

## Run simulations
set.seed(123)

nuisance_est_formula <- NA
nuisance_est_type <- "lm + optimal"
expect_y_est_type <- "lm"
e_est_type <- "lm"
expect_y_est_gam_withA <- NA

# 1. correct + correct -------------------------------------------------------
## ("Y ~ A*(X2 + t)") + ("R ~ X2 + t") -------------------------------------------------------
setting_list <- list(list(setting = 11, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ A*(X2 + t)")
e_formula <- as.formula("R ~ X2 + t")

for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
                             alpha0 = alpha0, lambda1 = lambda1,
                             alpha1 = alpha1, lambda2 = lambda2, functype = functype,
                             moderator = moderator,
                             expect_y_est_type = expect_y_est_type, expect_y_formula = expect_y_formula,
                             e_formula = e_formula, e_est_type = e_est_type,
                             expect_y_est_gam_withA = expect_y_est_gam_withA)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result, file_name)
  cat(file_name, "")
}

# 2. correct + mis -------------------------------------------------------
## ("Y ~  A*(X2 + t)") + ("R ~ t") -------------------------------------------------------
setting_list <- list(list(setting = 21, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~  A*(X2 + t)")
e_formula <- as.formula("R ~ t")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
                             alpha0 = alpha0, lambda1 = lambda1,
                             alpha1 = alpha1, lambda2 = lambda2, functype = functype,
                             moderator = moderator,
                             expect_y_est_type = expect_y_est_type, expect_y_formula = expect_y_formula,
                             e_formula = e_formula, e_est_type = e_est_type,
                             expect_y_est_gam_withA = expect_y_est_gam_withA)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result, file_name)
  cat(file_name, "")
}


# 3. mis + correct -------------------------------------------------------
## ("Y ~ X2 + t) + ("R ~ X2 + t) -------------------------------------------------------
setting_list <- list(list(setting = 31, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ X2 + t")
e_formula <- as.formula("R ~ X2 + t")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
                             alpha0 = alpha0, lambda1 = lambda1,
                             alpha1 = alpha1, lambda2 = lambda2, functype = functype,
                             moderator = moderator,
                             expect_y_est_type = expect_y_est_type, expect_y_formula = expect_y_formula,
                             e_formula = e_formula, e_est_type = e_est_type,
                             expect_y_est_gam_withA = expect_y_est_gam_withA)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result, file_name)
  cat(file_name, "")
}

# 4. mis + mis -------------------------------------------------------
## ("Y ~ X2 + t") + ("R ~ t") -------------------------------------------------------
setting_list <- list(list(setting = 41, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ X2 + t")
e_formula <- as.formula("R ~ t")

for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
                             alpha0 = alpha0, lambda1 = lambda1,
                             alpha1 = alpha1, lambda2 = lambda2, functype = functype,
                             moderator = moderator,
                             expect_y_est_type = expect_y_est_type, expect_y_formula = expect_y_formula,
                             e_formula = e_formula, e_est_type = e_est_type,
                             expect_y_est_gam_withA = expect_y_est_gam_withA)
  file_name <- paste0("sim", setting_list[[functype_i]]$setting, ".RDS")
  saveRDS(result, file_name)
  cat(file_name, "")
}

