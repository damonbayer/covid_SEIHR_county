# find parameters for overdispersion priors
# based on posteriors from NB spline
library(tidyverse)
library(rstan)
library(brms)
options(mc.cores = parallelly::availableCores())
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

labels <- c("hosp", "icu")

priors <- rbind(hosp_params, icu_params) %>%
          cbind(labels)

priors <- data.frame(priors)

colnames(priors) <- c("mean", "sd", "labels")
rownames(priors) <- NULL

file_name <- paste("overdisp_priors_countyid", indic, ".csv", sep = "")
write_csv(priors, paste("data/", file_name, sep = ""))
