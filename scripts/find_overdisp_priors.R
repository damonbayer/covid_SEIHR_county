library(tidyverse)
library(rstan)
library(brms)
library(fs)
library(cmdstanr)
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


# Cases -------------------------------------------------------------------
cases_fit <-
  brm(bf(est_cases ~ s(time)),
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

cases_fit_overdisp_draws <- as.vector(as_draws_array(cases_fit, variable = "shape"))
cases_overdisp_mean <- mean(log(cases_fit_overdisp_draws))
cases_overdisp_sd <- sd(log(cases_fit_overdisp_draws))

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


# Deaths ------------------------------------------------------------------
deaths_fit <-
  brm(bf(est_deaths ~ s(time)),
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

deaths_fit_overdisp_draws <- as.vector(as_draws_array(deaths_fit, variable = "shape"))
deaths_overdisp_mean <- mean(log(deaths_fit_overdisp_draws))
deaths_overdisp_sd <- sd(log(deaths_fit_overdisp_draws))

# Save Results ------------------------------------------------------------
priors <-
  tibble(
    datastream = c("cases", "hospitalized", "icu", "deaths"),
    mean = c(cases_overdisp_mean, hosp_overdisp_mean, icu_overdisp_mean, deaths_overdisp_mean),
    sd = c(cases_overdisp_sd, hosp_overdisp_sd, icu_overdisp_sd, deaths_overdisp_sd)
  )

dir_create("data/overdisp_priors")
file_name <- path("data/overdisp_priors", str_c("overdisp_priors_countyid=", target_county_id), ext = "csv")
write_csv(priors, file_name)
