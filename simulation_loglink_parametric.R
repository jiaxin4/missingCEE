library(rootSolve) 
library(mgcv)
library(tidyverse)
library(logbin)

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



beta_var_calculate_loglink <- function(dat, beta, moderator, 
                                       e_fit, expect_y_fit){
  # eta_hat for optimal formula
  S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
  beta <- as.matrix(beta, nrow = ncol(S_mat),ncol = 1)
  Sbeta <- as.numeric(S_mat %*% beta)
  
  dat$eta_hat <- (1-pt)*exp(-Sbeta)*dat$mu1_hat + pt*dat$mu0_hat
  i_total <- length(unique(dat$id))
  total_person_decisionpoint <- nrow(dat)
  
  #For e, logistic regression
  X_e <- model.matrix(e_fit) # Design matrix
  mu_e <- e_fit$fitted.values # Fitted values
  W_e <- diag(mu_e * (1 - mu_e))
  #For mu, linear regression
  dat <- dat %>% 
    mutate("AX2" = A*X2)
  variables_mu <- colnames(model.matrix(expect_y_fit))[-1] #instead of using model.matrix directly is becasue there are missing Ys.
  X_mu <- as.matrix(cbind(1, dat[, variables_mu]))
  W_mu <- diag(dat$expect_y_hat * (1 - dat$expect_y_hat))
  
  p <- ncol(X_e) + ncol(X_mu) + length(beta)
  meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
  for (it in 1:total_person_decisionpoint){
    
    #calculate derivative of Utilde wrt parameter of interest and nuisance parameter
    ##w.r.t beta
    dat_it <- dat[it, ]
    fx <- matrix(c(1, dat_it[, moderator]), nrow = length(beta), ncol = 1)
    dee_part1 <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$Y*(dat_it$A - pt)*dat_it$A)
    dee_part1 <- ifelse(is.na(dee_part1), 0, dee_part1)*fx%*%t(fx)
    dee_part2 <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$expect_y_hat*(dat_it$A - pt)*dat_it$A)*fx%*%t(fx)
    dee_Utilde_beta <- -(dat_it$R/dat_it$e_hat)*dee_part1 + ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*dee_part2
    ##w.r.t gamma_e
    gamma_e <- e_fit$coefficients
    dee_e_gamma_e <- as.numeric(exp(-gamma_e%*%X_e[it,])/(1 + exp(-gamma_e%*%X_e[it,]))^2)*(-X_e[it,])
    dee_Utilde_gamma_e <- (dat_it$R/dat_it$e_hat^2)*exp(-dat_it$A*t(fx)%*%beta)*(dat_it$Y - dat_it$expect_y_hat)*(dat_it$A - pt)
    dee_Utilde_gamma_e <- ifelse(is.na(dat_it$Y), 0, dee_Utilde_gamma_e)*fx%*%(dee_e_gamma_e)
    ##w.r.t gamma_mu
    gamma_mu <- expect_y_fit$coefficients
    dee_mu_gamma_mu <- as.numeric(exp(-gamma_mu%*%X_mu[it,])/(1 + exp(-gamma_mu%*%X_mu[it,]))^2)*(-X_mu[it,])
    dee_Utilde_gamma_mu <- as.numeric(((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*(dat_it$A - pt)*exp(-dat_it$A*t(fx)%*%beta))*fx%*%dee_mu_gamma_mu
    
    #calcualte derivative of nuisance parameter's own estimating equations
    ##dee of e
    dee_e <- -t(t(X_e[it,])) %*% W_e[it,it] %*% t(X_e[it,])
    ##dee of mu
    dee_mu <- -t(t(X_mu[it,])) %*% W_mu[it,it] %*% t(X_mu[it,])
    
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
    psi_gamma_mu <-  ifelse(is.na(dat_it$Y), 0, dat_it$Y - dat_it$expect_y_hat)*t(t(X_mu[it,]))
    
    ## EE of beta
    u1_it <- (exp(-dat_it$A*t(fx)%*%beta)*dat_it$Y - dat_it$eta_hat)*(dat_it$A - pt)
    u1_it <- ifelse(is.na(dat_it$Y), 0, as.numeric(u1_it))*fx
    u2_it <- as.numeric(exp(-dat_it$A*t(fx)%*%beta)*dat_it$expect_y_hat - dat_it$eta_hat)*(dat_it$A - pt)*fx
    ustar_it <- (dat_it$R/dat_it$e_hat)*u1_it - ((dat_it$R - dat_it$e_hat)/dat_it$e_hat)*u2_it
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
          #fit_gam_a0 <- gam(expect_y_formula, data = dat_a0)
          #fit_gam_a1 <- gam(expect_y_formula, data = dat_a1)
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
          #fit_gam <- gam(expect_y_formula, data = dat)
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
      beta_var <- beta_var_calculate_loglink(dat = dat, beta = beta_est_numerical, moderator = moderator,
                                             e_fit = e_fit, expect_y_fit = expect_y_fit)
      
      beta_est <- list(beta_est_numerical, beta_var)
      names(beta_est) <- c("beta_est_numerical", "asy_var")
      beta_est
    },
    ########################### If there is error (in fitting glm with loglink, then return NAs)
    error = function(cond) {
      message("\nCatched error")
      message(cond)
      S_mat <- as.matrix(cbind(rep(1, nrow(dat)), dat[, moderator]))
      beta_est <- list(beta_est_numerical = rep(NaN, ncol(S_mat)), 
                       asy_var = rep(NaN, ncol(S_mat)))
      names(beta_est) <- c("beta_est_numerical", "asy_var")
      return(beta_est)
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



### performance evaluation
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
    se_result <- apply(sqrt(asy_var_allsim), 1, mean)
    asy_var <- apply(asy_var_allsim, 1, mean)
    beta_numerical_result$se <- se_result
    
    
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
    
    # error rate
    nsim <- length(result_allsim)
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


#True parameter values
beta0 <- 1.1
beta1 <- -1.1
beta_true <- c(beta0, beta1)
alpha0 <- 0.3
lambda1 <- 0.2


t_total <- 20
pt = 0.4
nsim <- 1000

## Run simulations
set.seed(123)

i_total_list <- c(50, 100, 150, 200)


nuisance_est_type <- "lm + optimal"
nuisance_est_formula <- NA
expect_y_est_type <- "lm"
expect_y_est_gam_withA <- NA
e_est_type <- "lm"

# 1. correct + correct -------------------------------------------------------
setting_list <- list(list(setting = 11, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ A + X2 + AX2") #for creating interaction term in logbin
e_formula <- as.formula("R ~ X2")
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
setting_list <- list(list(setting = 21, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ A + X2 + AX2") #for creating interaction term in logbin
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
setting_list <- list(list(setting = 31, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ A") #for creating interaction term in logbin
e_formula <- as.formula("R ~ X2")
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
setting_list <- list(list(setting = 41, func_type = "linear", alpha1 = -0.5, lambda2 = 1.5))
moderator <- "X2"
expect_y_formula <- as.formula("Y ~ 1")
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



