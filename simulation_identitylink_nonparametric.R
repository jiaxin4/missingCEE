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




beta_var_calculate <- function(dat, beta, moderator){
  
  # eta_hat for optimal formula
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  Sbeta <- as.numeric(S_mat %*% beta)
  
  
  dat$eta_hat <- (pt + pt - 1)*Sbeta + (1-pt)*dat$mu1 + pt*dat$mu0
  
  i_total <- length(unique(dat$id))
  p <- length(beta)
  total_person_decisionpoint <- nrow(dat)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  for (it in 1:total_person_decisionpoint){
    dat_it <- dat[it, ]
    fx <- matrix(c(1, dat_it[, moderator]), nrow = p, ncol = 1)
    
    ## the derivative term
    dee_sum[it, , ] <- -((dat_it$A - pt)^2)*(fx%*%t(fx))           
    ## the meat term
    u1_it <- (dat_it$Y - dat_it$eta_hat - (dat_it$A - pt)*t(fx)%*%beta)*(dat_it$A - pt)
    u1_it <- ifelse(is.na(dat_it$Y), 0, as.numeric(u1_it))*fx
    u2_it <- (dat_it$A - pt)*fx*as.numeric(dat_it$expect_y_hat - dat_it$eta_hat - (dat_it$A - pt)*t(fx)%*%beta)
    ustar_it <- (dat_it$R/dat_it$e_hat)*u1_it - ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*u2_it
    meat_sum[it, ,] <- ustar_it%*%t(ustar_it)
  }
  dee <- apply(dee_sum, c(2,3), sum)/i_total
  dee_inv <- solve(dee)
  meat <- apply(meat_sum, c(2,3), sum) /i_total
  var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
  #var_cov <- dee_inv%*%meat%*%t(dee_inv)
  asy_var <- diag(var_cov)
  
  asy_var
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
    
    newdata_for_mu1_mu0 <-  data.frame(rbind(cbind(A = 1, moderator = dat[, moderator]),
                                             cbind(A = 0, moderator = dat[, moderator])))
    colnames(newdata_for_mu1_mu0) <- c("A", moderator)
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
  beta_var <- beta_var_calculate(dat = dat, beta = beta_est_numerical, moderator = moderator)
  
  beta_est <- list(beta_est_numerical, beta_var)
  names(beta_est) <- c("beta_est_numerical", "beta_var")
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





### Organize result
beta_result_function <- function(beta_hat_allsim, beta_true){
  beta_hat <- apply(beta_hat_allsim, 1, mean)
  bias <- apply(beta_hat_allsim - beta_true, 1, mean)
  mse <- apply((beta_hat_allsim - beta_true)^2, 1, mean)
  result <- list(beta_hat, bias, mse)
  names(result) <- c("beta_hat", "bias", "mse")
  result
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
      
      ci_in <- eval_ci(beta_hat$beta_est_numerical, beta_hat$beta_var, beta_true)
      
      #collect the result
      result <- list(beta_hat, ci_in, mean(dat$miss_prob_true))
      result_allsim[[isim]] <- result
    }
    beta_hat_allsim <-  (sapply(result_allsim, "[[", 1))[1,]
    beta_numerical_result <- beta_result_function(do.call(cbind, beta_hat_allsim), beta_true = beta_true)
    
    asy_var_allsim <- (sapply(result_allsim, "[[", 1))[2,]
    se_result <- apply(sqrt(do.call(cbind, asy_var_allsim)), 1, mean)
    beta_numerical_result$se_result <- se_result
    
    ci_in_allsim <- matrix(unlist(sapply(result_allsim, "[[", 2)), nrow = 3, byrow = FALSE)
    cp_result <- apply(ci_in_allsim, 1, mean)
    
    et_result <- mean(sapply(result_allsim, "[[", 3))
    
    hold <- list(i_total, t_total, 
                 beta_numerical_result, 
                 se_result,
                 cp_result,
                 et_result)
    names(hold) <- c("i_total", "t_total", 
                     "consistency_numerical", 
                     "se",
                     "coverage_probability", "missing_prob_true")
    result_i_total <- c(result_i_total, list(hold))
  }
  result_i_total
}


# Simulations
library(rootSolve) 
library(mgcv)
#True parameter values
beta0 <- 1.5
beta1 <- 2.1
beta_true <- c(beta0, beta1)
alpha0 <- 0.5
lambda1 <- 1.5
t_total <- 20
pt = 0.4
i_total_list <- c(50, 100, 150, 200)
nsim <- 1000


## Run simulations
set.seed(123)

# 1. correct + correct -------------------------------------------------------
## ("Y ~ s(X2) + s(t)", "byA") + ("R ~ s(X2) + s(t)") -------------------------------------------------------
setting_list <- list(list(setting = 11, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 12, func_type = "beta", alpha1 = -2, lambda2 = 1.5), 
                     list(setting = 13, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
nuisance_est_type <- "gam + optimal"
nuisance_est_formula <- NA
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2) + s(t)")
expect_y_est_gam_withA <- "byA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(X2) + s(t)")
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
## ("Y ~ s(X2) + s(t)", "byA) + ("R ~ s(t)") -------------------------------------------------------
setting_list <- list(list(setting = 21, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 22, func_type = "beta", alpha1 = -2, lambda2 = 1.5), 
                     list(setting = 23, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
nuisance_est_type <- "gam + optimal"
nuisance_est_formula <- NA
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2) + s(t)")
expect_y_est_gam_withA <- "byA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(t)")
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




# 4. mis + correct -------------------------------------------------------
## ("Y ~ s(t)", "byA") + ("R ~ s(X2) + s(t)") -------------------------------------------------------
setting_list <- list(list(setting = 61, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 62, func_type = "beta", alpha1 = -2, lambda2 = 1.5), 
                     list(setting = 63, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
nuisance_est_type <- "gam + optimal"
nuisance_est_formula <- NA
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(t)")
expect_y_est_gam_withA <- "byA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(X2) + s(t)")
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
## ("Y ~ s(X2) + s(t)", "notbyA") + ("R ~ s(t)") -------------------------------------------------------
setting_list <- list(list(setting = 41, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5), 
                     list(setting = 42, func_type = "beta", alpha1 = -2, lambda2 = 1.5), 
                     list(setting = 43, func_type = "periodic", alpha1 = 0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_est_type <- "gam"
expect_y_formula <- as.formula("Y ~ s(X2) + s(t)")
expect_y_est_gam_withA <- "notbyA"
e_est_type <- "gam"
e_formula <- as.formula("R ~ s(t)")

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