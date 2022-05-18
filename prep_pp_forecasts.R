library(tidyverse)
library(fs)

ci_widths <- c(0.5, 0.8, 0.95)
time_interval_in_days <- 7
county_id_key <- read_csv("data/county_id_key.csv")
county_region_key <- read_csv("data/county_region_key.csv")
county_id_region_key <- full_join(county_id_key, county_region_key) %>%
  rename(county_id = id)

time_date_key <-
  read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases) %>%
  distinct(time, date) %>%
  bind_rows(.,
            tibble(time = max(.$time) + 1:4,
                   date = max(.$date) + time_interval_in_days * (1:4)))

if (Sys.info()[["sysname"]] == "Linux") {
  results_dir <- "//dfs6/pub/bayerd/covid_SEIHR_county/results"
} else if (Sys.info()[["sysname"]] == "Darwin") {
  results_dir <- "results"
}

time_interval_in_days <- 7

posterior_predictive_samples <-
  path(results_dir, "posterior_predictive") %>%
  dir_ls() %>%
  enframe(name = NULL, value = "full_path") %>%
  mutate(file_name = full_path %>%
           path_file() %>%
           path_ext_remove()) %>%
  mutate(county_id = file_name %>%
           str_extract("(?<=county_id=)\\d+") %>%
           as.numeric()) %>%
  mutate(results = full_path %>%
           map(~read_csv(.) %>%
                 pivot_longer(-c(iteration, chain)) %>%
                 mutate(time = name %>%
                          str_extract("(?<=\\[)\\d+(?=\\])") %>%
                          as.numeric(),
                        name = name %>%
                          str_extract("^.+(?=\\[)") %>%
                          str_remove("data_")))) %>%
  select(county_id, results) %>%
  unnest(results)

posterior_predictive_summary_counties <-
  posterior_predictive_samples %>%
  select(-c(iteration, chain)) %>%
  group_by(county_id, name, time) %>%
  median_qi(.width = c(ci_widths, 0.9)) %>%
  left_join(county_id_region_key) %>%
  left_join(time_date_key) %>%
  select(county, date, name, everything(), -time, -county_id, -region) %>%
  arrange(county, date, name, .width)

posterior_predictive_summary_regions <-
  posterior_predictive_samples %>%
  left_join(county_id_region_key) %>%
  group_by(iteration, chain, name, time, region) %>%
  summarize(value = sum(value),
            .groups = "drop") %>%
  select(-c(iteration, chain)) %>%
  group_by(name, time, region) %>%
  median_qi(.width = ci_widths) %>%
  left_join(time_date_key) %>%
  select(region, date, name, everything(), -time) %>%
  arrange(region, date, name, .width)

posterior_predictive_summary_CA <-
  posterior_predictive_samples %>%
  select(-county_id) %>%
  group_by(iteration, chain, name, time) %>%
  summarize(value = sum(value),
            .groups = "drop") %>%
  select(-c(iteration, chain)) %>%
  group_by(name, time) %>%
  median_qi(.width = ci_widths) %>%
  left_join(time_date_key) %>%
  mutate(county = "California") %>%
  select(county, date, name, everything(), -time) %>%
  arrange(county, date, name, .width)

posterior_predictive_LEMMA_format <-
  bind_rows(posterior_predictive_summary_CA,
            posterior_predictive_summary_counties) %>%
  filter(.width == 0.9) %>%
  pivot_longer(cols = c(value, .lower, .upper),
               names_to = "value_type") %>%
  mutate(quantile =
           case_when(
             value_type == ".lower" ~ 0.05,
             value_type == "value" ~ 0.5,
             value_type == ".upper" ~ 0.95)) %>%
  select(county, date, quantile, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(cases = est_omicron_cases + est_other_cases) %>%
  select(date, county, quantile, hosp_census_with_covid = hospitalizations, cases)

write_csv(posterior_predictive_LEMMA_format, "posterior_predictive_LEMMA_format.csv")
write_csv(posterior_predictive_summary_CA, "results/posterior_predictive_summary_CA.csv")
write_csv(posterior_predictive_summary_regions, "results/posterior_predictive_summary_regions.csv")
write_csv(posterior_predictive_summary_counties %>% filter(.width != 0.9), "results/posterior_predictive_summary_counties.csv")
