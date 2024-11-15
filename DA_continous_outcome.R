library(tidyr)
library(dplyr)
library(MRTAnalysis)
library(knitr)
dat <- readRDS("suggestion_cleaned.RDS")

## check missing proportion
colSums(is.na(dat))

round(colSums(is.na(dat))/nrow(dat)*100, 2)

## comparison (1) complete case
dat_complete_case <- dat[which(!is.na(dat$logsteps30)),]

## comparison (2) handling missing data by imputing with 0
dat_imputezero <- dat
dat_imputezero[which(is.na(dat_imputezero$logsteps30)),"logsteps30"] <- 0 

## comparison (3) handling missing data by imputing with person-specific mean "
dat_imputemean <- dat %>%
  group_by(user) %>% 
  mutate(logsteps30 = ifelse(is.na(logsteps30), mean(logsteps30, na.rm = TRUE), logsteps30))
dat_imputemean <- ungroup(dat_imputemean)


# Estimand 1: fully marginal ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "logsteps30"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- c() #c() or "send"
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56)*send")
expectation_model_type <- "lm"
expectation_model_formula <- as.formula("logsteps30 ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56")
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## nonparametric nuisance models
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
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## comparators
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
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputezero,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56,     availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect


# Estimand 2: decision.index as moderator ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "logsteps30"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- "decision.index" #c() or "decision.index"
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log)*send")
expectation_model_type <- "lm"
expectation_model_formula <- as.formula("logsteps30 ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log") #already has A
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## nonparametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "logsteps30"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- "decision.index" #c() or "decision.index"
missingness_model_type <- "gam"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log, by = send) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
expectation_model_type <- "gam"
expectation_model_formula <- as.formula("logsteps30 ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator,
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## comparators
wcls_fit <- wcls(data = dat_complete_case,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputezero,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)

summary(wcls_fit)$causal_excursion_effect


# Estimand 3: weekday/weekend as moderator ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "logsteps30"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- "is.weekday"
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log)*send")
expectation_model_type <- "lm"
expectation_model_formula <- as.formula("logsteps30 ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log") #already has A
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## nonparametric nuisance models
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
result <- identitylink_cee_dr(dat = dat, 
                              id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator,
                              rand_prob = rand_prob, availability = availability, moderator = moderator,
                              missingness_model_type = missingness_model_type, 
                              missingness_model_formula = missingness_model_formula,
                              expectation_model_type = expectation_model_type,
                              expectation_model_formula = expectation_model_formula)
result

## comparators
wcls_fit <- wcls(data = dat_complete_case,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputezero,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "logsteps30",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)

summary(wcls_fit)$causal_excursion_effect














