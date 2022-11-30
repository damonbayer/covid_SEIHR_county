library(tidyverse)
source("src/plot_functions.R")

initialization_values <-
  read_csv("data/initialization_values.csv")

expit <- function(x)  1/(1+exp(-x))
logistic <- function(x)  1/(1+exp(-x))
logit <- function(x) -log(1/x - 1)

target_probabilities <- c(0.05, 0.5, 0.95)
n_samples <- 2000

remaining_prop_samples = sample(initialization_values$remaining_population / initialization_values$population, size = n_samples, replace = TRUE)

# E: 0.1% - 0.4%
E_prop_non_centered_mean <- logit(0.002)
E_prop_non_centered_sd <- 0.4
scales::percent(expit(qnorm(p = target_probabilities, mean = E_prop_non_centered_mean, sd = E_prop_non_centered_sd)))
E_prop_samples <- expit(rnorm(n = n_samples, mean = E_prop_non_centered_mean, sd = E_prop_non_centered_sd))

# I: 0.25% - 2%
I_prop_non_centered_mean <- logit(0.0075)
I_prop_non_centered_sd <- 0.6
scales::percent(expit(qnorm(p = target_probabilities, mean = I_prop_non_centered_mean, sd = I_prop_non_centered_sd)))
I_prop_samples <- expit(rnorm(n = n_samples, mean = I_prop_non_centered_mean, sd = I_prop_non_centered_sd))

# R: 5% - 20%
R_prop_non_centered_mean <- logit(0.1)
R_prop_non_centered_sd <- 0.5
scales::percent(expit(qnorm(p = target_probabilities, mean = R_prop_non_centered_mean, sd = R_prop_non_centered_sd)))
R_prop_samples <- expit(rnorm(n = n_samples, mean = R_prop_non_centered_mean, sd = R_prop_non_centered_sd))

# S: 80% - 95%
S_prop_samples <- remaining_prop_samples - (E_prop_samples + I_prop_samples + R_prop_samples)
scales::percent(quantile(S_prop_samples, target_probabilities))

# https://www.medrxiv.org/content/10.1101/2022.04.22.22274198v3
case_detection_rate_mean <- logit((1/6 + 1/8) / 2)
case_detection_rate_sd <- 0.1
rev(1/expit(qnorm(p = target_probabilities, mean = case_detection_rate_mean, sd = case_detection_rate_sd)))

# Previous paper
death_detection_rate_mean <- 2.3
death_detection_rate_sd <- 0.04
expit(qnorm(p = target_probabilities, mean = death_detection_rate_mean, sd = death_detection_rate_sd))

