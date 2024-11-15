# Code for paper "Doubly Robust Estimation of Causal Excursion Effects in Micro-Randomized Trials with Missing Longitudinal Outcomes"

Code to reproduce results in the paper [Doubly Robust Estimation of Causal Excursion Effects in Micro-Randomized Trials with Missing Longitudinal Outcomes] by Jiaxin Yu and Tianchen Qian.

Jiaxin Yu
2024.11.12

Section 1 provides a detailed guide on how to reproduce the simulation results in the paper and supplementary material. Section 2 provides how to conduct the data analysis in the paper.

## 1. Simulation Replication

The code for the simulations with identity link is with nonparametric estimation is [simulation_identitylink_nonparametric.R] (simulation_identitylink_nonparametric.R) (Section 5 in the paper). The code for the simulations with identity link is with parametric estimation and  [simulation_identitylink_parametric.R](simulation_identitylink_parametric.R). 

The code for the simulations with log link is [simulation_loglink_nonparametric.R](simulation_loglink_nonparametric.R) and  [simulation_loglink_parametric.R](simulation_loglink_parametric.R). 


## 2. Data analysis on HeartSteps

The data is publicly available at https://github.com/klasnja/HeartStepsV1/tree/main/data_files. The preprocessing code to prepare the dataset is [DA_preprocessing.R](DA_preprocessing.R).

#### Continuous outcome (Section 7 in the paper)

Code is [DA_continous_outcome.R](DA_continous_outcome.R) with functions from [DA_functions_continuous_outcome.R](DA_functions_continuous_outcome.R). 

### Binary outcome (Section E in the supplementary material)

Code is [DA_binary_outcome.R](DA_binary_outcome.R) with functions from [DA_functions_binary_outcome.R](DA_functions_binary_outcome.R). 








