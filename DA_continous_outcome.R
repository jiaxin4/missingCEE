rm(list = ls())
library(MRTAnalysis)
library(knitr)
dat <- readRDS("suggestion_cleaned.RDS")
source("DA_functions_continuous_outcome.R")
library(rootSolve) 
library(mgcv)
library(tidyverse)
expit <- function(x){
  1/(1+exp(-x))
}


params_list <- list(list("original_missing"), 
                    list(int = -2, send = 0.3, dp = 0.008, weekday = -0.6),
                    list(int = -2, send = 1.3, dp = 0.008, weekday = -0.6),
                    list(int = -2, send = 2.3, dp = 0.008, weekday = -0.6)
)

###########################

miss_percent_list <- c(0, 10, 20, 30)

dat_marginal <- dat_moderated_time <- dat_moderated_weekday <- c()
for (i_miss in 1:length(params_list)){
    
  if (i_miss == 1){
    dat_new <- dat
  } else {
    int <- params_list[[i_miss]]$int
    send <- params_list[[i_miss]]$send
    dp <- params_list[[i_miss]]$dp
    weekday <- params_list[[i_miss]]$weekday
    ## adding missingness
    prob_missing <- expit(int + send*dat$send + dp*dat$decision.index + weekday*dat$is.weekday)
    artificial_missing = rbinom(length(dat$steps30), size = 1, prob = prob_missing)
    dat$steps30_new = ifelse(artificial_missing == 1 | is.na(dat$steps30), NA, dat$steps30)
    
    ### look at the distribution of missingness
    #summary(dat$steps30_new)
    cat(length(which(is.na(dat$steps30_new)))/nrow(dat)*100, "% is currently missing")
    dat %>% 
      group_by(user) %>% 
      summarise(missing_steps = sum(is.na(steps30_new))) %>% 
      ggplot(aes(x = missing_steps)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
      labs(title = "Histogram of Missing Observations per User",
           x = "Number of Missing Observations",
           y = "Count of Users") +
      theme_bw()
    
    dat_new <- dat %>% 
      dplyr::select(user, 
                    decision.index = decision.index,
                    avail = avail,
                    send = send,
                    steps30 = steps30_new,
                    is.weekday = is.weekday,
                    homework.location = homework.location,
                    steps30pre.zero = steps30pre.zero) %>% 
      mutate(logsteps30 = log(steps30 + 0.5),
             steps30_observe_indicator = ifelse(is.na(steps30), 0, 1),
             steps30pre.zero = ifelse(steps30_observe_indicator, steps30pre.zero, 0),
             steps30pre.log = log(steps30pre.zero + 0.5),
             steps30pre.zero.bin_56 = ifelse(steps30pre.zero >= 56, 1, 0))
  }
  
  ## comparison (1) complete case
  dat_complete_case <- dat_new[which(!is.na(dat_new$logsteps30)),]
  ## comparison (2) handling missing data by imputing with 0
  dat_imputezero <- dat_new
  dat_imputezero[which(is.na(dat_imputezero$logsteps30)),"logsteps30"] <- 0 
  ## comparison (3) handling missing data by imputing with person-specific mean "
  dat_imputemean <- dat_new %>%
    group_by(user) %>% 
    mutate(logsteps30 = ifelse(is.na(logsteps30), mean(logsteps30, na.rm = TRUE), logsteps30))
  dat_imputemean <- ungroup(dat_imputemean)
  
  ## Marginal -----------------------------------------------------------------------
  ## Nuisance models: nonparametric
  id <- "user"
  treatment <- "send"
  outcome <- "logsteps30"
  response_observe_indicator <- "steps30_observe_indicator"
  rand_prob <- 0.6
  availability <- "avail"
  moderator <- c() #c() or "send"
  missingness_model_type <- "gam"
  missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log, by = send) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
  expectation_model_type <- "gam"
  expectation_model_formula <- as.formula("logsteps30 ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
  result <- identitylink_cee_dr(dat = dat_new, 
                                id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                                rand_prob = rand_prob, availability = availability, moderator = moderator,
                                missingness_model_type = missingness_model_type, 
                                missingness_model_formula = missingness_model_formula,
                                expectation_model_type = expectation_model_type,
                                expectation_model_formula = expectation_model_formula)
  beta_hat <- result$beta_hat
  ci_low <- result$ci$ci_low
  ci_high <- result$ci$ci_high
  dr_beta_colored_nonparam <- cbind(beta_hat, ci_low, ci_high)
  ## comparison (1) complete case
  wcls_fit <- wcls(data = dat_complete_case,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ 1,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,
                   availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  cc_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (2) handling missing data by imputing with 0
  wcls_fit <- wcls(data = dat_imputezero,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ 1,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputezero_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (3) handling missing data by imputing with person-specific mean "
  wcls_fit <- wcls(data = dat_imputemean,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ 1,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputemean_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  tab <- data.frame(rbind(dr_beta_colored_nonparam,
                          cc_beta_colored, imputezero_beta_colored, imputemean_beta_colored))
  tab$method <- rep(c("Doubly_Robust_nonparam", "Complete_Case", "Imputed_zero", "Imputed_mean"), each = length(hold[,1]))
  tab$beta <- rep(paste0("beta", 0:(length(hold[,1]) - 1)), times =4)
  
  miss_percent <- round(length(which(is.na(dat$steps30_new)))/nrow(dat)*100, 1)
  tab$miss_percent <- miss_percent
  tab$miss_percent_tick <- miss_percent_list[i_miss]
  
  dat_marginal <- rbind(dat_marginal, tab)
  
  
  ## moderated by time -----------------------------------------------------------------------
  ## Nuisance models: nonparametric
  id <- "user"
  treatment <- "send"
  outcome <- "logsteps30"
  response_observe_indicator <- "steps30_observe_indicator"
  rand_prob <- 0.6
  availability <- "avail"
  moderator <- "decision.index"
  missingness_model_type <- "gam"
  missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log, by = send) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
  expectation_model_type <- "gam"
  expectation_model_formula <- as.formula("logsteps30 ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
  result <- identitylink_cee_dr(dat = dat_new, 
                                id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                                rand_prob = rand_prob, availability = availability, moderator = moderator,
                                missingness_model_type = missingness_model_type, 
                                missingness_model_formula = missingness_model_formula,
                                expectation_model_type = expectation_model_type,
                                expectation_model_formula = expectation_model_formula)
  beta_hat <- result$beta_hat
  ci_low <- result$ci$ci_low
  ci_high <- result$ci$ci_high
  dr_beta_colored_nonparam <- cbind(beta_hat, ci_low, ci_high)
  ## comparison (1) complete case
  wcls_fit <- wcls(data = dat_complete_case,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ decision.index,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,
                   availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  cc_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (2) handling missing data by imputing with 0
  wcls_fit <- wcls(data = dat_imputezero,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ decision.index,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputezero_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (3) handling missing data by imputing with person-specific mean "
  wcls_fit <- wcls(data = dat_imputemean,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ decision.index,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputemean_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  tab <- data.frame(rbind(dr_beta_colored_nonparam,
                          cc_beta_colored, imputezero_beta_colored, imputemean_beta_colored))
  tab$method <- rep(c("Doubly_Robust_nonparam", "Complete_Case", "Imputed_zero", "Imputed_mean"), each = length(hold[,1]))
  tab$beta <- rep(paste0("beta", 0:(length(hold[,1]) - 1)), times =4)
  
  miss_percent <- round(length(which(is.na(dat$steps30_new)))/nrow(dat)*100, 1)
  tab$miss_percent <- miss_percent
  tab$miss_percent_tick <- miss_percent_list[i_miss]
  
  dat_moderated_time <- rbind(dat_moderated_time, tab)
  
  
  ## moderated by weekday -----------------------------------------------------------------------
  ## Nuisance models: nonparametric
  id <- "user"
  treatment <- "send"
  outcome <- "logsteps30"
  response_observe_indicator <- "steps30_observe_indicator"
  rand_prob <- 0.6
  availability <- "avail"
  moderator <- "is.weekday"
  missingness_model_type <- "gam"
  missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log, by = send) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
  expectation_model_type <- "gam"
  expectation_model_formula <- as.formula("logsteps30 ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
  result <- identitylink_cee_dr(dat = dat_new, 
                                id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                                rand_prob = rand_prob, availability = availability, moderator = moderator,
                                missingness_model_type = missingness_model_type, 
                                missingness_model_formula = missingness_model_formula,
                                expectation_model_type = expectation_model_type,
                                expectation_model_formula = expectation_model_formula)
  beta_hat <- result$beta_hat
  ci_low <- result$ci$ci_low
  ci_high <- result$ci$ci_high
  dr_beta_colored_nonparam <- cbind(beta_hat, ci_low, ci_high)
  ## comparison (1) complete case
  wcls_fit <- wcls(data = dat_complete_case,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ is.weekday,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,
                   availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  cc_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (2) handling missing data by imputing with 0
  wcls_fit <- wcls(data = dat_imputezero,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ is.weekday,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputezero_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  ## comparison (3) handling missing data by imputing with person-specific mean "
  wcls_fit <- wcls(data = dat_imputemean,
                   id = "user",
                   outcome = "logsteps30",
                   treatment = "send",
                   rand_prob = 0.6,
                   moderator_formula = ~ is.weekday,
                   control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                   numerator_prob = 0.6
  )
  hold <- summary(wcls_fit)$causal_excursion_effect
  imputemean_beta_colored <- cbind(hold[,1], hold[,2], hold[,3])
  
  tab <- data.frame(rbind(dr_beta_colored_nonparam,
                          cc_beta_colored, imputezero_beta_colored, imputemean_beta_colored))
  tab$method <- rep(c("Doubly_Robust_nonparam", "Complete_Case", "Imputed_zero", "Imputed_mean"), each = length(hold[,1]))
  tab$beta <- rep(paste0("beta", 0:(length(hold[,1]) - 1)), times =4)
  
  miss_percent <- round(length(which(is.na(dat$steps30_new)))/nrow(dat)*100, 1)
  tab$miss_percent <- miss_percent
  tab$miss_percent_tick <- miss_percent_list[i_miss]
  
  dat_moderated_weekday <- rbind(dat_moderated_weekday, tab)
  
}


# marginal figure ------------------------------------------------------
method_color <- c("Doubly_Robust_nonparam" = "#009E73",
                  "Complete_Case" = "#CC79A7",
                  "Imputed_mean" = "#56B4E9",
                  "Imputed_zero" = "#E69F00"
)

dat_marginal$method <- factor(dat_marginal$method, levels = c("Doubly_Robust_nonparam", 
                                                              "Complete_Case",
                                                              "Imputed_mean",
                                                              "Imputed_zero"))
miss_xaxis_labels <- paste(sort(unique(dat_marginal$miss_percent)), "%", sep = "")


custom_theme <- theme(axis.text.x = element_text(size = 12),
                      axis.title.x = element_text(size = 13),
                      axis.text.y = element_text(size = 12),
                      axis.title.y = element_text(size = 13))

p_marginal <- dat_marginal %>% 
  filter(beta == "beta0",
         method != "Doubly_Robust_param") %>% 
  ggplot(aes(x = miss_percent_tick, y = beta_hat, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = method),
                linewidth = 1, 
                width = 5, 
                position = position_dodge(width = 5)) +
  theme_bw() +
  scale_x_continuous(
    breaks = miss_percent_list,       # Replace these with your actual numeric values
    labels = miss_xaxis_labels
  ) +
  xlab("Proportion of missing observations") +
  ylab("Intercept") +
  custom_theme +
  scale_color_manual(name = "Methods", 
                     values = method_color,
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero")) +
  scale_shape_manual(name = "Methods", 
                     values = c(16, 17, 15, 18),
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero"))

# moderated by time figure ------------------------------------------------------
dat_moderated_time$method <- factor(dat_moderated_time$method, levels = c("Doubly_Robust_nonparam", 
                                                                          "Complete_Case",
                                                                          "Imputed_mean",
                                                                          "Imputed_zero"))
miss_xaxis_labels <- paste(sort(unique(dat_moderated_time$miss_percent)), "%", sep = "")

p_moderated_time_beta0 <- dat_moderated_time %>% 
  filter(beta == "beta0",
         method != "Doubly_Robust_param") %>% 
  ggplot(aes(x = miss_percent_tick, y = beta_hat, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = method),
                linewidth = 1, 
                width = 5, 
                position = position_dodge(width = 5)) +
  theme_bw() +
  scale_x_continuous(
    breaks = miss_percent_list,       # Replace these with your actual numeric values
    labels = miss_xaxis_labels
  ) +
  xlab("Proportion of missing observations") +
  ylab("Intercept") +
  custom_theme +
  scale_color_manual(name = "Methods", 
                     values = method_color,
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero")) +
  scale_shape_manual(name = "Methods", 
                     values = c(16, 17, 15, 18),
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero"))

p_moderated_time_beta1 <- dat_moderated_time %>% 
  filter(beta == "beta1",
         method != "Doubly_Robust_param") %>% 
  ggplot(aes(x = miss_percent_tick, y = beta_hat, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = method),
                linewidth = 1, 
                width = 5, 
                position = position_dodge(width = 5)) +
  theme_bw() +
  scale_x_continuous(
    breaks = miss_percent_list,      
    labels = miss_xaxis_labels
  ) +
  xlab("Proportion of missing observations") +
  ylab("Slope (decision point)") +
  custom_theme +
  scale_color_manual(name = "Methods", 
                     values = method_color,
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero")) +
  scale_shape_manual(name = "Methods", 
                     values = c(16, 17, 15, 18),
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero"))

# moderated by weekday figure ------------------------------------------------------
dat_moderated_weekday$method <- factor(dat_moderated_weekday$method, levels = c("Doubly_Robust_nonparam", 
                                                                          "Complete_Case",
                                                                          "Imputed_mean",
                                                                          "Imputed_zero"))
miss_xaxis_labels <- paste(sort(unique(dat_moderated_weekday$miss_percent)), "%", sep = "")

p_moderated_weekday_beta0 <- dat_moderated_weekday %>% 
  filter(beta == "beta0",
         method != "Doubly_Robust_param") %>% 
  ggplot(aes(x = miss_percent_tick, y = beta_hat, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = method),
                linewidth = 1, 
                width = 5, 
                position = position_dodge(width = 5)) +
  theme_bw() +
  scale_x_continuous(
    breaks = miss_percent_list,       # Replace these with your actual numeric values
    labels = miss_xaxis_labels
  ) +
  xlab("Proportion of missing observations") +
  ylab("Intercept") +
  custom_theme +
  scale_color_manual(name = "Methods", 
                     values = method_color,
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero")) +
  scale_shape_manual(name = "Methods", 
                     values = c(16, 17, 15, 18),
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero"))

p_moderated_weekday_beta1 <- dat_moderated_weekday %>% 
  filter(beta == "beta1",
         method != "Doubly_Robust_param") %>% 
  ggplot(aes(x = miss_percent_tick, y = beta_hat, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = method),
                linewidth = 1, 
                width = 5, 
                position = position_dodge(width = 5)) +
  theme_bw() +
  scale_x_continuous(
    breaks = miss_percent_list,      
    labels = miss_xaxis_labels
  ) +
  xlab("Proportion of missing observations") +
  ylab("Slope (is.weekday)") +
  custom_theme  +
  scale_color_manual(name = "Methods", 
                     values = method_color,
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero")) +
  scale_shape_manual(name = "Methods", 
                     values = c(16, 17, 15, 18),
                     labels = c("Doubly_Robust_nonparam" =  "DR-nonparametric",
                                "Complete_Case" = "Complete-case", 
                                "Imputed_mean" = "Imputed-mean", "Imputed_zero" = "Imputed-zero"))

#combine them -------------------------------------------------------------------
library(cowplot)
# Extract legend from one plot
shared_legend <- get_legend(
  p_marginal + 
    theme(legend.position = "right",
                     legend.text = element_text(size = 12),
                     legend.title = element_text(size = 12)) 
)

# Remove legends from all individual plots
p_marginal_clean <- p_marginal + theme(legend.position = "none") +  labs(x = NULL) 
p_moderated_time_beta0_clean <- p_moderated_time_beta0 + theme(legend.position = "none") +  labs(x = NULL) 
p_moderated_time_beta1_clean <- p_moderated_time_beta1 + theme(legend.position = "none") +  labs(x = NULL) 
p_moderated_weekday_beta0_clean <- p_moderated_weekday_beta0 + theme(legend.position = "none") 
p_moderated_weekday_beta1_clean <- p_moderated_weekday_beta1 + theme(legend.position = "none")

# Arrange the plots: First row with marginal plot and legend, second row with two moderated plots
first_column <- plot_grid(p_marginal_clean,
                          p_moderated_time_beta0_clean,
                          p_moderated_weekday_beta0_clean,
                          nrow = 3, ncol = 1,
                          rel_widths = c(1,1,1),
                          rel_heights = c(1,1,1),
                          align = "v",
                          labels = c("(a)", "(b)", "(c)"),
                          label_size = 18)
second_column <- plot_grid(shared_legend,
                          plot_grid(p_moderated_time_beta1_clean,
                          p_moderated_weekday_beta1_clean,
                          nrow = 2, ncol = 1,
                          rel_widths = c(1,1),
                          rel_heights = c(1,1),
                          align = "v"),
                          ncol = 1, nrow = 2,
                          rel_heights = c(1, 2))
p <- plot_grid(first_column, second_column, 
          nrow = 1, ncol = 2,
          rel_heights = c(1,1))

p














