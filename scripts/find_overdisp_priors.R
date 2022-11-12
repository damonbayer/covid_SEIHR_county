# find parameters for overdispersion priors
# based on posteriors from NB spline
library(tidyverse)
library(rstan)
library(brms)
library(fs)
options(mc.cores = parallelly::availableCores(),
        brms.backend = "cmdstanr")
rstan_options(auto_write = TRUE)

source("src/spline_functions.R")

# command args for array job ----------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
indic <- as.integer(args[1])

county_id_key <- read_csv("data/county_id_key.csv")
county_id <- county_id_key %>% filter(id == indic) %>% pull(county)

data <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
        filter(county == county_id)

# find hospitalization priors ---------------------------------------------
hosp_spline <- run_nb_spline(data = data,
                             response = "hosp")
hosp_overdisp <- choose_kappa_params(hosp_spline)
hosp_params <- hosp_overdisp$par

# find icu priors ---------------------------------------------------------
icu_spline <- run_nb_spline(data, response = "icu")
icu_overdisp <- choose_kappa_params(icu_spline)
icu_params <- icu_overdisp$par

# save the results --------------------------------------------------------
priors <-
  tibble(param_type = c("mean", "sd"), hosp = hosp_params, icu = icu_params) %>%
  pivot_longer(-param_type, names_to = "datastream") %>%
  pivot_wider(names_from = param_type, values_from = value)

file_name <- path("data/overdisp_priors", str_c("overdisp_priors_countyid=", indic), ext = "csv")
write_csv(priors, file_name)
