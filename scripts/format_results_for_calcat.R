library(tidyverse)
library(fs)
results_dir <- "results"

initialization_values <- read_csv("data/initialization_values.csv")

tidy_posterior_predictive <- read_csv(path(results_dir, "tidy_posterior_predictive", ext = "csv"))
tidy_vector_posterior_generated_quantities <- read_csv(path(results_dir, "tidy_vector_posterior_generated_quantities", ext = "csv"))

generated_quantities_LEMMA_format <-
  tidy_vector_posterior_generated_quantities %>%
  filter(name == "Rₜ_t") %>%
  select(-c(.point, .interval)) %>%
  filter(.width == 0.90) %>%
  pivot_longer(cols = c(value, .lower, .upper), names_to = "value_type") %>%
  mutate(quantile = case_when(
    value_type == "value" ~ 0.5,
    value_type == ".lower" ~ (1 - .width) / 2,
    value_type == ".upper" ~ 1 - (1 - .width) / 2
  )) %>%
  mutate(quantile = round(quantile, digits = 2)) %>%
  select(-c(.width, value_type)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename(Rt = "Rₜ_t" )

posterior_predictive_LEMMA_format <-
  tidy_posterior_predictive %>%
  select(-c(.point, .interval)) %>%
  filter(.width == 0.90) %>%
  pivot_longer(cols = c(value, .lower, .upper), names_to = "value_type") %>%
  mutate(quantile = case_when(
    value_type == "value" ~ 0.5,
    value_type == ".lower" ~ (1 - .width) / 2,
    value_type == ".upper" ~ 1 - (1 - .width) / 2
  )) %>%
  select(-c(.width, value_type)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename(hosp_census_with_covid = hospitalizations,
         deaths = new_deaths,
         cases = new_cases,
         cumulative_new_deaths = cumulative_deaths) %>%
  left_join(initialization_values %>%
              select(county, initial_cumulative_deaths = D)) %>%
  mutate(cumulative_deaths = cumulative_new_deaths + initial_cumulative_deaths) %>%
  select(-c(cumulative_new_deaths, initial_cumulative_deaths))

results_calcat_format <-
  left_join(posterior_predictive_LEMMA_format,
            generated_quantities_LEMMA_format) %>%
  left_join(initialization_values %>%
              select(county, id)) %>%
  arrange(id) %>%
  select(-id)

write_csv(results_calcat_format, "results_calcat_format.csv")
