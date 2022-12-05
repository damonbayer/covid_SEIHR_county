# https://twitter.com/trvrb/status/1599830501803843584
library(tidyverse)
library(ggdist)
expit <- function(x)  1/(1+exp(-x))
logistic <- function(x)  1/(1+exp(-x))
logit <- function(x) -log(1/x - 1)

case_detection_rate_non_centered_mean = -2.9
case_detection_rate_non_centered_sd = 0.3

1 / expit(qnorm(p = c(0.05, 0.5, 0.95), mean = case_detection_rate_non_centered_mean, case_detection_rate_non_centered_sd))


# I -> H -> ICU -> D
# should be 0.04%-0.07%

IHR_non_centered_mean = -4.3
IHR_non_centered_sd = 0.25

HICUR_non_centered_mean = -1.69
HICUR_non_centered_sd = 0.2

ICUDR_non_centered_mean = -1.59
ICUDR_non_centered_sd = 0.2


n_samples <- 20000

IHR_samples <- expit(rnorm(n_samples, IHR_non_centered_mean, IHR_non_centered_sd))

IHHICUR_samples <- expit(rnorm(n_samples, HICUR_non_centered_mean, HICUR_non_centered_sd))

ICUDR_samples <- expit(rnorm(n_samples, ICUDR_non_centered_mean, ICUDR_non_centered_sd))

# Current values seem about right
tibble(IFR = IHR_samples * IHHICUR_samples * ICUDR_samples) %>%
  ggplot(aes(IFR)) +
  stat_halfeye() +
  scale_x_continuous(labels = scales::percent)
