# function for data analysis
library(logbin)
expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(p){
  log(p/(1-p))
}

loglink_cee_dr <- function(dat, 
                          id, treatment, outcome, response_observe_indicator, rand_prob, availability, moderator,
                          missingness_model_type, missingness_model_formula,
                          expectation_model_type, expectation_model_formula){
  i_total <- length(unique(dat[[id]]))
  Y <- dat[[outcome]]
  A <- dat[[treatment]]
  R <- dat[[response_observe_indicator]]
  pt <- rand_prob
  I <- dat[[availability]]
  S_mat <- as.matrix( cbind( rep(1, nrow(dat)), dat[, moderator]))

  #Step 1: Fit the nuisance models
  
  ## Fit the missingness model P(Rt = 1|Ht, At)
  if (missingness_model_type == "glm"){
    #missingness_model_formula <- as.formula(paste0(response_observe_indicator, " ~ ", paste0(missingness_model_vars, "", collapse = " + ")))
    miss_fit <- glm(missingness_model_formula,  data = dat, family = binomial(link = "logit"))
    r_hat <- predict(miss_fit, newdata = dat, type = "response")
  } else if (missingness_model_type == "gam"){
    #missingness_model_formula <- as.formula(paste0(response_observe_indicator, " ~ ", paste0("s(", missingness_model_vars, ")", collapse = " + ")))
    miss_fit <- gam(missingness_model_formula, data = dat, family = binomial(link = "logit"))
    r_hat <- predict(miss_fit, newdata = dat, type = "response")
  }

  ## Fit the expectation model E(Yt,1|Ht, At, It = 1)
  if (expectation_model_type == "glm"){
    dat_forfit <- cbind(Y, A, dat[,attr(terms(expectation_model_formula), "term.labels")],
                        dat[,attr(terms(expectation_model_formula), "term.labels")]*A)
    colnames(dat_forfit) <- c(outcome, treatment, attr(terms(expectation_model_formula), "term.labels"),
                              paste(treatment, attr(terms(expectation_model_formula), "term.labels"), sep = "_"))
    logbin_formula <- paste0(outcome, "~", paste(colnames(dat_forfit)[-1], collapse = "+"))
    initial_values <- coef(logbin(logbin_formula, data = dat_forfit))

    expectation_model_formula_withA <- paste0(outcome, " ~ ", paste("send*", attr(terms(expectation_model_formula), "term.labels"), collapse = " + "))
    expect_y_fit <- glm(expectation_model_formula_withA, family=binomial(link="log"),
                        start = rep(0.1, ncol(dat_forfit)),
                        data = dat, weights = dat[[availability]])
    expect_y_hat <- predict(expect_y_fit, newdata = dat, type = "response")

    newdata_for_mu1_mu0 <- dat[, attr(terms(expectation_model_formula), "term.labels")]
    newdata_for_mu1_mu0 <- rbind(cbind(A = 1, newdata_for_mu1_mu0),
                                 cbind(A = 0, newdata_for_mu1_mu0))
    colnames(newdata_for_mu1_mu0) <- c(treatment, attr(terms(expectation_model_formula), "term.labels"))
    mu_s <- cbind(newdata_for_mu1_mu0,
                  expect_y = predict(expect_y_fit, newdata = newdata_for_mu1_mu0, type = "response"))
    mu1_hat <- mu_s[mu_s$A == 1, "expect_y"]
    mu0_hat <- mu_s[mu_s$A == 0, "expect_y"]

  } else if (expectation_model_type == "gam"){
    dat_a0 <- subset(dat, A == 0)
    dat_a1 <- subset(dat, A == 1)
    dat_a0$I_a0 <- dat_a0[[availability]]
    dat_a1$I_a1 <- dat_a1[[availability]]
    fit_gam_a0 <- gam(expectation_model_formula, family=binomial(link="log"), data = dat_a0, weights = I_a0)
    fit_gam_a1 <- gam(expectation_model_formula, family=binomial(link="log"), data = dat_a1, weights = I_a1)

    dat_set_all_a_to_0 <- dat
    dat_set_all_a_to_0$A <- 0
    dat_set_all_a_to_0$I_a0 <- dat_set_all_a_to_0[[availability]]
    mu0_hat <- predict(fit_gam_a0, newdata = dat_set_all_a_to_0, type = "response")

    dat_set_all_a_to_1 <- dat
    dat_set_all_a_to_1$A <- 1
    dat_set_all_a_to_1$I_a1 <- dat_set_all_a_to_1[[availability]]
    mu1_hat <- predict(fit_gam_a1, newdata = dat_set_all_a_to_1, type = "response")

    expect_y_hat <- ifelse(A, mu1_hat, mu0_hat)
  }
  # if (missingness_model_type == "glm" & expectation_model_type == "glm"){
  #   nuisancefit_glm <- readRDS("nuisancefit_glm.RDS")
  #   miss_fit <- nuisancefit_glm[[1]][[1]]
  #   r_hat <- nuisancefit_glm[[1]][[2]]
  #   expect_y_fit <- nuisancefit_glm[[2]][[1]]
  #   expect_y_hat <- nuisancefit_glm[[2]][[2]]
  #   mu1_hat <- nuisancefit_glm[[2]][[3]]
  #   mu0_hat <- nuisancefit_glm[[2]][[4]]
  # } else if (missingness_model_type == "gam" & expectation_model_type == "gam"){
  #   nuisancefit_glm <- readRDS("nuisancefit_gam.RDS")
  #   miss_fit <- nuisancefit_glm[[1]][[1]]
  #   r_hat <- nuisancefit_glm[[1]][[2]]
  #   expect_y_hat <- nuisancefit_glm[[2]][[1]]
  #   mu1_hat <- nuisancefit_glm[[2]][[2]]
  #   mu0_hat <- nuisancefit_glm[[2]][[3]]
  # }
  
  
  #Step 2: Solve the estimating equation
  ## Estimate beta
  beta_numerical_solution_loglink <- function(){
    
    ee <- function(beta){
      Sbeta <- as.numeric(S_mat %*% beta)
      
      eta_hat <- (1-pt)*exp(-Sbeta)*mu1_hat + pt*mu0_hat
      
      U1 <- (exp(-A*Sbeta)*Y - eta_hat)*(A - pt)*I
      U1 <- ifelse(is.na(U1), 0, U1)
      EU <- (exp(-A*Sbeta)*expect_y_hat - eta_hat)*I
      
      Ustar <- rep(NA, length(beta)) 
      for (ibeta in 1:length(beta)) {
        Ustar[ibeta] <- sum((R/r_hat)*U1 * S_mat[, ibeta] 
                            - ((R - r_hat)/r_hat)*(A - pt)*S_mat[, ibeta]*EU)
      }
      
      Ustar <- Ustar/i_total
      return(Ustar)
    }
    solution <- multiroot(ee, rep(0, ncol(S_mat)),  useFortran = FALSE)
    beta_hat_solve <- solution$root
    
    beta_hat_solve
  } 

  beta_hat <- beta_numerical_solution_loglink()
  
  ## Estimate the variance
  beta_var_calculate_lm_loglink <- function(beta = beta_hat){
    
    # eta_hat for optimal formula
    Sbeta <- as.numeric(S_mat %*% beta)
    eta_hat <- (1-pt)*exp(-Sbeta)*mu1_hat + pt*mu0_hat
    total_person_decisionpoint <- length(Y)
    
    #For missingness model, logistic regression
    X_e <- model.matrix(miss_fit) # Design matrix
    W_e <- diag(r_hat * (1 - r_hat))
    #For expectation model, logistic regression
    X_mu <- as.matrix(cbind(1, A, dat[,attr(terms(expectation_model_formula), "term.labels")],
                  dat[,attr(terms(expectation_model_formula), "term.labels")]*A))
    W_mu <- diag(expect_y_hat * (1 - expect_y_hat))

    p <- ncol(X_e) + ncol(X_mu) + length(beta)
    meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
    dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
    for (it in 1:total_person_decisionpoint){
      
      #calculate derivative of Utilde wrt parameter of interest and nuisance parameter
      ## w.r.t beta
      fx <- matrix(S_mat[it,])
      dee_part1 <- as.numeric(exp(-A[it]*t(fx)%*%beta)*Y[it]*(A[it] - pt)*A[it])
      dee_part1 <- ifelse(is.na(dee_part1), 0, dee_part1)*fx%*%t(fx)
      dee_part2 <- as.numeric(exp(-A[it]*t(fx)%*%beta)*expect_y_hat[it]*(A[it] - pt)*A[it])*fx%*%t(fx)
      dee_Utilde_beta <- -(R[it]/r_hat[it])*dee_part1*I[it] + ((R[it] - r_hat[it])/r_hat[it])*dee_part2*I[it]
      ## w.r.t gamma_e
      gamma_e <- miss_fit$coefficients
      dee_e_gamma_e <- as.numeric(exp(-gamma_e%*%X_e[it,])/(1 + exp(-gamma_e%*%X_e[it,]))^2)*(-X_e[it,])
      dee_Utilde_gamma_e <- (R[it]/r_hat[it]^2)*exp(-A[it]*t(fx)%*%beta)*(Y[it] - expect_y_hat[it])*(A[it] - pt)*I[it]
      dee_Utilde_gamma_e <- ifelse(is.na(Y[it]), 0, dee_Utilde_gamma_e)*fx%*%(dee_e_gamma_e)
      
      ## w.r.t gamma_mu
      gamma_mu <- expect_y_fit$coefficients
      dee_mu_gamma_mu <- as.numeric(exp(-gamma_mu%*%X_mu[it,])/(1 + exp(-gamma_mu%*%X_mu[it,]))^2)*(-X_mu[it,])
      dee_Utilde_gamma_mu <- as.numeric(((R[it] - r_hat[it])/r_hat[it])*(A[it] - pt)*exp(-A[it]*t(fx)%*%beta))*I[it]*fx%*%dee_mu_gamma_mu
      
      #calcualte derivative of nuisance parameter's own estimating equations
      ## dee of e
      dee_e <- -t(t(X_e[it,])) %*% W_e[it,it] %*% t(X_e[it,])
      ## dee of mu
      dee_mu <- -t(t(X_mu[it,])) %*% W_mu[it,it] %*% t(X_mu[it,]) * I[it]
      
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
      psi_gamma_e <-  (R[it] - r_hat[it])*t(t(X_e[it,]))
      ## EE of gamma_mu
      psi_gamma_mu <-  ifelse(is.na(Y[it]), 0, Y[it] - expect_y_hat[it])*t(t(X_mu[it,]))*I[it]

      ## EE of beta
      u1_it <- (exp(-A[it]*t(fx)%*%beta)*Y[it] - eta_hat[it])*(A[it] - pt)
      u1_it <- ifelse(is.na(Y[it]), 0, as.numeric(u1_it))*fx*I[it]
      u2_it <- as.numeric(exp(-A[it]*t(fx)%*%beta)*expect_y_hat[it] - eta_hat[it])*(A[it] - pt)*fx*I[it]
      ustar_it <- (R[it]/r_hat[it])*u1_it - ((R[it] - r_hat[it])/r_hat[it])*u2_it
      meat_sum[it, ,] <- rbind(psi_gamma_e, psi_gamma_mu, ustar_it)%*%t(rbind(psi_gamma_e, psi_gamma_mu, ustar_it))
      
    }
    dee <- apply(dee_sum, c(2,3), sum)/i_total
    dee_inv <- solve(dee)
    meat <- apply(meat_sum, c(2,3), sum) /i_total
    var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
    asy_var <- diag(var_cov)
    asy_var_beta <- asy_var[(ncol(X_e) + ncol(X_mu) + 1):p]
    
    asy_var_beta
    
  }
  #####################
  beta_var_calculate_gam_loglink <- function(beta = beta_hat){
    # eta_hat for optimal formula
    Sbeta <- as.numeric(S_mat %*% beta)
    eta_hat <- (1-pt)*exp(-Sbeta)*mu1_hat + pt*mu0_hat
    total_person_decisionpoint <- length(Y)
    
    p <- length(beta)
    meat_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
    dee_sum <- array(NA, dim = c(total_person_decisionpoint, p, p))
    for (it in 1:total_person_decisionpoint){
      fx <- matrix(S_mat[it,])    
      
      ## the derivative term
      dee_part1 <- as.numeric(exp(-A[it]*t(fx)%*%beta)*Y[it]*(A[it] - pt)*A[it])
      dee_part1 <- ifelse(is.na(dee_part1), 0, dee_part1)*fx%*%t(fx)*I[it]
      dee_part2 <- as.numeric(exp(-A[it]*t(fx)%*%beta)*expect_y_hat[it]*(A[it] - pt)*A[it])*fx%*%t(fx)*I[it]
      dee_sum[it, , ] <- -(R[it]/r_hat[it])*dee_part1 + ((R[it] - r_hat[it])/r_hat[it])*dee_part2     
      
      ## the meat term
      u1_it <- (exp(-A[it]*t(fx)%*%beta)*Y[it] - eta_hat[it])*(A[it] - pt)
      u1_it <- ifelse(is.na(Y[it]), 0, as.numeric(u1_it))*fx*I[it]
      u2_it <- as.numeric(exp(-A[it]*t(fx)%*%beta)*expect_y_hat[it] - eta_hat[it])*(A[it] - pt)*fx*I[it]
      ustar_it <- (R[it]/r_hat[it])*u1_it - ((R[it] - r_hat[it])/r_hat[it])*u2_it
      meat_sum[it, ,] <- ustar_it%*%t(ustar_it)
    }
    dee <- apply(dee_sum, c(2,3), sum)/i_total
    dee_inv <- solve(dee)
    meat <- apply(meat_sum, c(2,3), sum) /i_total
    var_cov <- dee_inv%*%meat%*%t(dee_inv)/i_total
    asy_var <- diag(var_cov)
    
    asy_var
  }
  
  
  #################
  if (missingness_model_type == "glm" & expectation_model_type == "glm"){
    beta_var <- beta_var_calculate_lm_loglink()
  } else if (missingness_model_type == "gam" & expectation_model_type == "gam"){
    beta_var <- beta_var_calculate_gam_loglink()
  }
  
  ci_low <- beta_hat - 1.96*sqrt(beta_var)
  ci_high <- beta_hat + 1.96*sqrt(beta_var)
  ci <- list(ci_low, ci_high)
  names(ci) <- c("ci_low", "ci_high")
  
  beta_est <- list(beta_hat, beta_var, ci)
  names(beta_est) <- c("beta_hat", "asy_var", "ci")
  beta_est
}

















