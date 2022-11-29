library(tidyverse)
library(rstan)
library(brms)
library(fs)
options(
  mc.cores = parallelly::availableCores(),
  brms.backend = "cmdstanr"
)
rstan_options(auto_write = TRUE)

# Read Data ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
target_county_id <- as.integer(args[1])


county_id_pop <- read_csv("data/county_id_pop.csv")
target_county_name <- county_id_pop %>%
  filter(id == target_county_id) %>%
  pull(county)

data <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  filter(county == target_county_name)

# Model Arguments ---------------------------------------------------------
seed <- 17
cores <- parallelly::availableCores()
iter <- 10000
warmup <- 5000
thin <- 10
refresh <- 0
adapt_delta <- 0.99

# Hospitalization ---------------------------------------------------------
hosp_fit <-
  brm(bf(hospitalized ~ s(time)),
    data = data,
    family = negbinomial(),
    cores = cores,
    seed = seed,
    iter = iter,
    warmup = warmup,
    thin = thin,
    refresh = refresh,
    control = list(adapt_delta = adapt_delta)
  )

hosp_fit_overdisp_draws <- as.vector(as_draws_array(hosp_fit, variable = "shape"))
hosp_overdisp_mean <- mean(log(hosp_fit_overdisp_draws))
hosp_overdisp_sd <- sd(log(hosp_fit_overdisp_draws))


# ICU ---------------------------------------------------------------------
icu_fit <-
  brm(bf(icu ~ s(time)),
    data = data,
    family = negbinomial(),
    cores = cores,
    seed = seed,
    iter = iter,
    warmup = warmup,
    thin = thin,
    refresh = refresh,
    control = list(adapt_delta = adapt_delta)
  )

icu_fit_overdisp_draws <- as.vector(as_draws_array(icu_fit, variable = "shape"))
icu_overdisp_mean <- mean(log(icu_fit_overdisp_draws))
icu_overdisp_sd <- sd(log(icu_fit_overdisp_draws))

# Save Results ------------------------------------------------------------
priors <-
  tibble(
    datastream = c("hospitalized", "icu"),
    mean = c(hosp_overdisp_mean, icu_overdisp_mean),
    sd = c(hosp_overdisp_sd, icu_overdisp_sd)
  )

dir_create("data/overdisp_priors")
file_name <- path("data/overdisp_priors", str_c("overdisp_priors_countyid=", target_county_id), ext = "csv")
write_csv(priors, file_name)
