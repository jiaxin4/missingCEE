# Functions 
even_indicator <- function(num){
  ifelse(num %% 2 == 0, 1, 0)
}

gen_nonlinear <- function(whichtype, xt, t, t_total, alpha0, lambda1){
  if (whichtype == "beta"){
    out <- alpha0 + lambda1*(dbeta(xt/6, 2, 2)) 
  } else if (whichtype == "periodic"){
    out <- alpha0 + lambda1*sin(xt)
  } else if (whichtype == "linear"){
    out <- alpha0 + lambda1*xt
  }
  out
}

## Data generating functions
dgm_nonlinear_loglink <- function(i_total, t_total, beta0, beta1, 
                                  alpha0, lambda1, alpha1, lambda2, 
                                  functype){
  dat <- c()
  for (i in 1:i_total){
    #Generate the data:
    t <- 1:t_total
    #X1i <- runif(t_total, 0, 2) #time-varying covariate
    X2i <- runif(t_total, 0, 2)
    Ai <- rbinom(t_total, 1, pt) #treatment variable
    
    delta_i <- beta0 + beta1*X2i#CEE
    #delta_i <- beta0 
    eta_i <- gen_nonlinear(whichtype = functype, xt = X2i, t = t, t_total = t_total, 
                           alpha0 = alpha0, lambda1 = lambda1)
    miss_prob_genmod <- gen_nonlinear(whichtype = functype, xt = X2i, t = t, t_total = t_total,
                                      alpha0 = alpha1, lambda1 = lambda2)
    
    pYi <- eta_i*exp(Ai*delta_i) 
    Yi <- rbinom(t_total, 1, pYi) #binary outcome
    
    #Missing data:
    miss_prob_true <- exp(miss_prob_genmod)/(1 + exp(miss_prob_genmod))
    Ri <- 1*rbinom(t_total, 1, miss_prob_true)
    Yi <- ifelse(Ri == 1, Yi, NA)
    
    #So the data here is {Ai, X2i, RiYi, Yi}
    dat <- rbind(dat, cbind(i, t, Ai, X2i, Ri, Yi, miss_prob_true, pYi))
  }
  colnames(dat) <- c("id", "t", "A", "X2", "R", "Y", "miss_prob_true", "pY")
  dat <- as.data.frame(dat)
  dat
}




## Estimate beta
############################################################# numerical solution
beta_numerical_solution_loglink <- function(dat, moderator){
  i_total <- length(unique(dat$id))
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  
  ee <- function(beta){
    Sbeta <- as.numeric(S_mat %*% beta)
    
    dat$eta_hat <- (1-pt)*exp(-Sbeta)*dat$mu1_hat + pt*dat$mu0_hat
    
    U1 <- (exp(-(dat$A)*Sbeta)*dat$Y - dat$eta_hat)*(dat$A - pt)
    U1 <- ifelse(is.na(U1), 0, U1)
    EU <- exp(-(dat$A)*Sbeta)*dat$expect_y_hat - dat$eta_hat 
    
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





beta_var_calculate_loglink <- function(dat, beta, moderator){
  
  # eta_hat for optimal formula
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  beta <- as.matrix(beta, nrow = ncol(S_mat))
  Sbeta <- as.numeric(S_mat %*%beta)
  
  dat$eta_hat <- (1-pt)*exp(-Sbeta)*dat$mu1_hat + pt*dat$mu0_hat
  
  i_total <- length(unique(dat$id))
  p <- length(beta)
  total_person_decisionpoint <- nrow(dat)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  for (it in 1:total_person_decisionpoint){
    dat_it <- dat[it, ]
    fx <- as.matrix(rbind(rep(1, nrow(dat_it)), dat_it[, moderator]))
    
    ## the derivative term
    dee_part1 <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$Y*(dat_it$A - pt)*dat_it$A)
    dee_part1 <- ifelse(is.na(dee_part1), 0, dee_part1)*fx%*%t(fx)
    dee_part2 <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$expect_y_hat*(dat_it$A - pt)*dat_it$A)*fx%*%t(fx)
    dee_sum[it, , ] <- -(dat_it$R/dat_it$e_hat)*dee_part1 + ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*dee_part2     
    
    ## the meat term
    u1_it <- (exp(-dat_it$A*t(fx)%*%beta)*dat_it$Y - dat_it$eta_hat)*(dat_it$A - pt)
    u1_it <- ifelse(is.na(dat_it$Y), 0, as.numeric(u1_it))*fx
    u2_it <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$expect_y_hat - dat_it$eta_hat)*(dat_it$A - pt)*fx
    ustar_it <- (dat_it$R/dat_it$e_hat)*u1_it - ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*u2_it
    meat_sum[it, ,] <- ustar_it%*%t(ustar_it)
  }
  dee <- solve(apply(dee_sum, c(2,3), sum)/i_total)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee%*%meat%*%t(dee)/i_total
  asy_var <- diag(var_cov)
  
  asy_var
}

## Estimate beta
beta_est_function_loglink <- function(dat, moderator,
                                      e_formula, expect_y_est_type, expect_y_formula, 
                                      ...){
  i_total <- length(unique(dat$id))
  
  ##############
  beta_solution <- tryCatch(
    ######################### If there is no error, then return est_beta
    {
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
      # browser()
      if (expect_y_est_type == "lm"){
        
        dat$AX2 <- dat$A*dat$X2
        initial_values <- coef(logbin(expect_y_formula, data = dat))
        expect_y_fit <- glm(expect_y_formula, family=binomial(link="log"), data = dat,
                            start = initial_values)
        expect_y_hat <- predict(expect_y_fit, newdata = dat, type = "response")
        dat$expect_y_hat <- expect_y_hat
        
        newdata_for_mu1_mu0 <-  data.frame(rbind(cbind(A = 1, t = dat$t, moderator = dat$X2),
                                                 cbind(A = 0, t = dat$t, moderator = dat$X2)))
        colnames(newdata_for_mu1_mu0) <- c("A", "t", "X2")
        newdata_for_mu1_mu0$AX2 <- newdata_for_mu1_mu0$A*newdata_for_mu1_mu0$X2
        
        mu_s <- cbind(newdata_for_mu1_mu0, expect_y = predict(expect_y_fit, newdata = newdata_for_mu1_mu0, type = "response"))
        mu1_hat <- mu_s[mu_s$A == 1, "expect_y"]
        mu0_hat <- mu_s[mu_s$A == 0, "expect_y"]
        dat$mu1_hat <- mu1_hat
        dat$mu0_hat <- mu0_hat
      } else if (expect_y_est_type == "gam"){
        if (expect_y_est_gam_withA == "byA"){
          dat_a0 <- subset(dat, A == 0)
          dat_a1 <- subset(dat, A == 1)
          fit_gam_a0 <- gam(expect_y_formula, family=binomial(link="log"), data = dat_a0)
          fit_gam_a1 <- gam(expect_y_formula, family=binomial(link="log"), data = dat_a1)
          
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
          fit_gam <- gam(expect_y_formula, family=binomial(link="log"), data = dat)
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
      beta_est_numerical <- beta_numerical_solution_loglink(dat = dat, moderator = moderator)
      
      ## Estimate the variance
      beta_var <- beta_var_calculate_loglink(dat = dat, beta = beta_est_numerical, moderator = moderator)
      
      beta_est <- list(beta_est_numerical, beta_var)
      names(beta_est) <- c("beta_est_numerical", "asy_var")
      beta_est
    },
    ########################### If there is error (in fitting glm with loglink, then return NAs)
    error = function(cond) {
      message("\nCatched error")
      message(cond)
      S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
      return(list(beta_est_numerical = rep(NaN, ncol(S_mat)), 
                  asy_var = rep(NaN, ncol(S_mat))))
    }
  )
  
  
}



eval_ci <- function(beta_hat, beta_var, beta_true){
  ci_low <- beta_hat - 1.96*sqrt(beta_var)
  ci_high <- beta_hat + 1.96*sqrt(beta_var)
  result <- list()
  for (beta_i in 1:length(beta_hat)){
    beta_i_in <- ifelse(beta_true[beta_i] <= ci_high[beta_i] & beta_true[beta_i] >= ci_low[beta_i], 1, 0)
    result[[beta_i]] <- beta_i_in
  }
  # beta_0_in <- ifelse(beta_true[1] <= ci_high[1] & beta_true[1] >= ci_low[1], 1, 0)
  # beta_1_in <- ifelse(beta_true[2] <= ci_high[2] & beta_true[2] >= ci_low[2], 1, 0)
  beta_in <- ifelse(all(unlist(result) == 1), 1, 0)
  result <- c(result, beta_in)
  result
}



# simulation function
sim_nonlineardgm_loglink <- function(i_total_list, t_total, beta0, beta1, 
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
      dat <- dgm_nonlinear_loglink(i_total = i_total, t_total = t_total,  
                                   beta0 = beta0, beta1 = beta1, 
                                   alpha0 = alpha0, lambda1 = lambda1,
                                   alpha1 = alpha1, lambda2 = lambda2,
                                   functype = functype)
      
      beta_hat <- beta_est_function_loglink(dat = dat, moderator = moderator,
                                            e_formula = e_formula, 
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


### Organize result
beta_result_function <- function(beta_hat_allsim, beta_true){
  beta_hat <- apply(beta_hat_allsim, 1, mean, na.rm=TRUE)
  bias <- apply(beta_hat_allsim - beta_true, 1, mean, na.rm=TRUE)
  mse <- apply((beta_hat_allsim - beta_true)^2, 1, mean, na.rm=TRUE)
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
    se_result <- apply(sqrt(asy_var_allsim), 1, mean, na.rm = TRUE)
    asy_var <- apply(asy_var_allsim, 1, mean, na.rm = TRUE)
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
    cp_result <- apply(ci_in_allsim, 1, mean, na.rm = TRUE)
    # CP based on true variance
    true_varaince <-  apply(beta_hat_allsim, 1, var, na.rm = TRUE)
    ci_truevariance <- list()
    for (i in 1:ncol(beta_hat_allsim)){
      ci_truevariance[[i]] <- eval_ci(beta_hat = beta_hat_allsim[,i],
                                      beta_var = true_varaince, 
                                      beta_true = beta_true)
    }
    cp_result_truevariance <- c(mean(sapply(ci_truevariance, "[[", 1), na.rm = TRUE),
                                mean(sapply(ci_truevariance, "[[", 2), na.rm = TRUE))
    
    beta_numerical_result$asy_var <- asy_var
    beta_numerical_result$true_var <- true_varaince
    
    #True observed response probability
    et_result <- mean(sapply(result_allsim, "[[", 3))
    
    error_prop <- mean(rowSums(is.na(beta_hat_allsim))/nsim)
    
    #collect results
    hold <- list(i_total, t_total,
                 beta_numerical_result,
                 cp_result,
                 et_result,
                 cp_result_truevariance,
                 error_prop)
    names(hold) <- c("i_total", "t_total",
                     "beta_numerical_result",
                     "coverage_probability", 
                     "missing_prob_true",
                     "cp_result_truevariance",
                     "error_prop")
    result_organized <- c(result_organized, list(hold))
  }
  result_organized
}





## Simulation
library(rootSolve) 
library(mgcv)
library(tidyverse)
library(logbin)

#True parameter values
beta0 <- 1.1
beta1 <- -1.1
beta_true <- c(beta0, beta1)
alpha0 <- 0.3
lambda1 <- 0.2


t_total <- 20
pt = 0.4
i_total_list <- c(50, 100, 150, 200)
nsim <- 1000

## Run simulations
set.seed(123)


# 1. correct + correct -------------------------------------------------------
## ("Y ~ s(X2) by A") + ("R ~ s(X2)") -------------------------------------------------------
setting_list <- list(list(setting = 11, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 12, func_type = "beta", alpha1 = -0.2, lambda2 = 1.5), 
                     list(setting = 13, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2)") 
expect_y_est_gam_withA <- "byA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(X2)")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm_loglink(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
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
## ("Y ~ s(X2) by A") + ("R ~ 1") -------------------------------------------------------
setting_list <- list(list(setting = 21, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 22, func_type = "beta", alpha1 = -0.2, lambda2 = 1.5), 
                     list(setting = 23, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2)") 
expect_y_est_gam_withA <- "byA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ 1")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm_loglink(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
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
## ("Y ~ s(X2) notbyA") + ("R ~ s(X2)") -------------------------------------------------------
setting_list <- list(list(setting = 31, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 32, func_type = "beta", alpha1 = -0.2, lambda2 = 1.5), 
                     list(setting = 33, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2)") 
expect_y_est_gam_withA <- "notbyA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(X2)")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm_loglink(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
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
## ("Y ~ s(X2)") + ("R ~ 1") -------------------------------------------------------
setting_list <- list(list(setting = 41, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 42, func_type = "beta", alpha1 = -0.2, lambda2 = 1.5), 
                     list(setting = 43, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2)") 
expect_y_est_gam_withA <- "notbyA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ 1")
for (functype_i in 1:length(setting_list)){
  functype <- setting_list[[functype_i]]$func_type
  alpha1 <- setting_list[[functype_i]]$alpha1
  lambda2 <- setting_list[[functype_i]]$lambda2
  
  result <- sim_nonlineardgm_loglink(i_total_list = i_total_list, t_total = t_total, beta0 = beta0, beta1 = beta1,
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

