library(tidyverse)
library(tidybayes)
library(fs)
library(furrr)
# Setup -------------------------------------------------------------------
print(str_c(parallelly::availableCores(), " cores detected"))
plan(strategy = multisession, workers = parallelly::availableCores())
results_dir <- "results"

ci_widths <- c(0.5, 0.7, 0.9)
time_interval_in_days <- 7
n_forecast_times <- 12
county_id_pop <- read_csv("data/county_id_pop.csv")
county_region_key <- read_csv("data/county_region_key.csv")

dat <- read_csv("data/cases_hospitalizations_by_county.csv")

time_date_key <-
  dat %>%
  distinct(time, date = end_date) %>%
  add_row(., time = 0, date = min(.[["date"]]) - time_interval_in_days, .before = 1) %>%
  add_row(.,
    time = max(.[["time"]]) + 1:n_forecast_times,
    date = max(.[["date"]]) + time_interval_in_days * 1:n_forecast_times
  )


# Tidy functions ----------------------------------------------------------
tidy_generated_quantities_file <- function(file_path) {
  read_csv(file_path) %>%
    select(-c(iteration, chain)) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    median_qi(.width = ci_widths) %>%
    mutate(
      time = name %>%
        str_extract("(?<=\\[)\\d+(?=\\])") %>%
        as.numeric()
    ) %>%
    mutate(name = if_else(
      str_detect(name, "\\[\\d+\\]"),
      name %>%
        str_extract("^.+(?=\\[)") %>%
        str_remove("data_"),
      name
    )) %>%
    mutate(time = if_else(
      str_ends(name, "_mean", negate = T),
      time - 1,
      time
    )) %>%
    left_join(time_date_key) %>%
    select(name, date, everything(), -time) %>%
    arrange(name, date)
}

tidy_generated_quantities_dir <- function(dir_path) {
  all_generated_quantities <-
    dir_ls(dir_path) %>%
    enframe(name = NULL, value = "full_path") %>%
    mutate(id = full_path %>%
      str_extract("(?<=county_id=)\\d+") %>%
      as.numeric()) %>%
    arrange(id) %>%
    left_join(county_id_pop %>% select(-population)) %>%
    mutate(results = future_map(full_path, tidy_generated_quantities_file)) %>%
    select(-full_path, -id) %>%
    unnest(results)

  scalar_generated_quantities <-
    all_generated_quantities %>%
    filter(is.na(date)) %>%
    select(-date)

  vector_generated_quantities <-
    all_generated_quantities %>%
    filter(!is.na(date))

  list(
    scalar_generated_quantities = scalar_generated_quantities,
    vector_generated_quantities = vector_generated_quantities
  )
}

tidy_predictive_file <- function(file_path) {
  all_predictive_wide <- read_csv(file_path)

  cumulative_deaths_summary <-
    all_predictive_wide %>%
    select(iteration, chain, starts_with("data_new_deaths")) %>%
    pivot_longer(-c(iteration, chain)) %>%
    mutate(
      time = name %>%
        str_extract("(?<=\\[)\\d+(?=\\])") %>%
        as.numeric()
    ) %>%
    mutate(name = if_else(
      str_detect(name, "\\[\\d+\\]"),
      name %>%
        str_extract("^.+(?=\\[)") %>%
        str_remove("data_"),
      name
    )) %>%
    group_by(iteration, chain) %>%
    mutate(value = cumsum(value)) %>%
    ungroup() %>%
    mutate(name = "cumulative_deaths") %>%
    select(-c(iteration, chain)) %>%
    group_by(name, time) %>%
    median_qi(.width = ci_widths)

  non_cumulative_predictive_summary <-
    all_predictive_wide %>%
    select(-c(iteration, chain)) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    median_qi(.width = ci_widths) %>%
    mutate(
      time = name %>%
        str_extract("(?<=\\[)\\d+(?=\\])") %>%
        as.numeric()
    ) %>%
    mutate(name = if_else(
      str_detect(name, "\\[\\d+\\]"),
      name %>%
        str_extract("^.+(?=\\[)") %>%
        str_remove("data_"),
      name
    ))


  bind_rows(
    non_cumulative_predictive_summary,
    cumulative_deaths_summary
  ) %>%
    left_join(time_date_key) %>%
    select(name, date, everything(), -time) %>%
    arrange(name, date)
}

tidy_predictive_dir <- function(dir_path) {
  all_generated_quantities <-
    dir_ls(dir_path) %>%
    enframe(name = NULL, value = "full_path") %>%
    mutate(id = full_path %>%
      str_extract("(?<=county_id=)\\d+") %>%
      as.numeric()) %>%
    arrange(id) %>%
    left_join(county_id_pop %>% select(-population)) %>%
    mutate(results = future_map(full_path, tidy_predictive_file)) %>%
    select(-full_path, -id) %>%
    unnest(results)
}


# Tidy directories --------------------------------------------------------
tidy_prior_generated_quantities <- tidy_generated_quantities_dir(path(results_dir, "prior_generated_quantities"))
tidy_posterior_generated_quantities <- tidy_generated_quantities_dir(path(results_dir, "posterior_generated_quantities"))

tidy_prior_predictive <- tidy_predictive_dir(path(results_dir, "prior_predictive"))
tidy_posterior_predictive <- tidy_predictive_dir(path(results_dir, "posterior_predictive"))

# Save files --------------------------------------------------------------
write_csv(tidy_prior_generated_quantities$scalar_generated_quantities, path(results_dir, "tidy_scalar_prior_generated_quantities", ext = "csv"))
write_csv(tidy_prior_generated_quantities$vector_generated_quantities, path(results_dir, "tidy_vector_prior_generated_quantities", ext = "csv"))

write_csv(tidy_posterior_generated_quantities$scalar_generated_quantities, path(results_dir, "tidy_scalar_posterior_generated_quantities", ext = "csv"))
write_csv(tidy_posterior_generated_quantities$vector_generated_quantities, path(results_dir, "tidy_vector_posterior_generated_quantities", ext = "csv"))

write_csv(tidy_prior_predictive, path(results_dir, "tidy_prior_predictive", ext = "csv"))
write_csv(tidy_posterior_predictive, path(results_dir, "tidy_posterior_predictive", ext = "csv"))
