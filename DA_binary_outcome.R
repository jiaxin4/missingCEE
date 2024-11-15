library(MRTAnalysis)
library(knitr)
dat <- readRDS("suggestion_cleaned.RDS")

## dichotomize by mean
dat$steps30_binary_median <- ifelse(dat$steps30 > median(dat$steps30, na.rm = TRUE), 1, 0)

## comparison (1) complete case
dat_complete_case <- dat[which(!is.na(dat$logsteps30)),]

## comparison (2) handling missing data by imputing with person-specific mean "
dat_imputemean <- dat %>%
  group_by(user) %>% 
  mutate(steps30_binary_median = ifelse(is.na(steps30_binary_median), round(mean(steps30_binary_median, na.rm = TRUE)), steps30_binary_median))
dat_imputemean <- ungroup(dat_imputemean)

## comparison (3) impute all missing as 0
dat_imputezero <- dat
dat_imputezero[is.na(dat_imputezero$steps30_binary_median), "steps30_binary_median"] <- 0

## comparison (4) impute all missing as 1
dat_imputeone <- dat
dat_imputeone[is.na(dat_imputeone$steps30_binary_median), "steps30_binary_median"] <- 1


# Estimand 1: fully marginal ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- c() #c() or "send"
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56)*send")
expectation_model_type <- "glm"
expectation_model_formula <- as.formula("steps30_binary_median ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log")
result <- loglink_cee_dr(dat = dat, 
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
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- c() #c() or "send"
missingness_model_type <- "gam"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
expectation_model_type <- "gam"
expectation_model_formula <- as.formula("steps30_binary_median ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
result <- loglink_cee_dr(dat = dat, 
                         id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                         rand_prob = rand_prob, availability = availability, moderator = moderator,
                         missingness_model_type = missingness_model_type, 
                         missingness_model_formula = missingness_model_formula,
                         expectation_model_type = expectation_model_type,
                         expectation_model_formula = expectation_model_formula)


## comparators
wcls_fit <- wcls(data = dat_complete_case,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputezero,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputeone,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ 1,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect


# Estimand 2: decision index as moderator ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- "decision.index" 
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56)*send")
expectation_model_type <- "glm"
expectation_model_formula <- as.formula("steps30_binary_median ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log")
result <- loglink_cee_dr(dat = dat, 
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
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- "decision.index" 
missingness_model_type <- "gam"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
expectation_model_type <- "gam"
expectation_model_formula <- as.formula("steps30_binary_median ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
result <- loglink_cee_dr(dat = dat, 
                         id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                         rand_prob = rand_prob, availability = availability, moderator = moderator,
                         missingness_model_type = missingness_model_type, 
                         missingness_model_formula = missingness_model_formula,
                         expectation_model_type = expectation_model_type,
                         expectation_model_formula = expectation_model_formula)



## comparators
wcls_fit <- wcls(data = dat_complete_case,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "steps30_binary_median",
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
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputeone,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~ decision.index,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect


# Estimand 3: decision index as moderator ----------------------------------------------
## parametric nuisance models
id <- "user"
treatment <- "send"
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- c("is.weekday")
missingness_model_type <- "glm"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ (decision.index + homework.location + is.weekday + steps30pre.log + steps30pre.zero.bin_56)*send")
expectation_model_type <- "glm"
expectation_model_formula <- as.formula("steps30_binary_median ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log")
result <- loglink_cee_dr(dat = dat, 
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
outcome <- "steps30_binary_median"
response_observe_indicator <- "steps30_observe_indicator"
rand_prob <- 0.6
availability <- "avail"
moderator <- c("is.weekday")
missingness_model_type <- "gam"
missingness_model_formula <- as.formula("steps30_observe_indicator ~ s(decision.index, by = send) + s(steps30pre.log) + homework.location*send + is.weekday*send + steps30pre.zero.bin_56*send")
expectation_model_type <- "gam"
expectation_model_formula <- as.formula("steps30_binary_median ~ s(decision.index) + homework.location + is.weekday + steps30pre.zero.bin_56 + s(steps30pre.log)")
result <- loglink_cee_dr(dat = dat, 
                         id = id, treatment = treatment, outcome = outcome, response_observe_indicator = response_observe_indicator, 
                         rand_prob = rand_prob, availability = availability, moderator = moderator,
                         missingness_model_type = missingness_model_type, 
                         missingness_model_formula = missingness_model_formula,
                         expectation_model_type = expectation_model_type,
                         expectation_model_formula = expectation_model_formula)



## comparators
wcls_fit <- wcls(data = dat_complete_case,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~  is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputemean,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~  is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputezero,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~  is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect

wcls_fit <- wcls(data = dat_imputeone,
                 id = "user",
                 outcome = "steps30_binary_median",
                 treatment = "send",
                 rand_prob = 0.6,
                 moderator_formula = ~  is.weekday,
                 control_formula = ~ decision.index + homework.location + is.weekday + steps30pre.zero.bin_56 + steps30pre.log,
                 availability = "avail",
                 numerator_prob = 0.6
)
summary(wcls_fit)$causal_excursion_effect



