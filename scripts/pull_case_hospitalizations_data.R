library(tidyverse)
library(ckanr)
library(lubridate)
library(splines)
library(fs)

time_interval_in_days <- 7

results_dir <- "results"

case_reporting_delay_ecdf <- read_rds("src/case_reporting_delay_ecdf.rds")
death_reporting_delay_ecdf <- read_rds("src/death_delay_ecdf.rds")
county_region_key <- read_csv("data/county_region_key.csv")

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

ckanr_setup(url = "https://data.ca.gov")
ckan <- quiet(ckanr::src_ckan("https://data.ca.gov"))

# get resources
resources <- rbind(
  resource_search("name:covid-19", as = "table")$results,
  resource_search("name:hospitals by county", as = "table")$results
)


cases_deaths_url <- resources %>%
  filter(name == "Statewide COVID-19 Cases Deaths Tests") %>%
  pull(url)

hosp_url <- resources %>%
  filter(name == "Statewide Covid-19 Hospital County Data") %>%
  pull(url)

cases_raw <-
  read_csv(cases_deaths_url) %>%
  filter(area_type == "County") %>%
  rename(county = area) %>%
  filter(!(county %in% c("Out of state", "Unknown")))

county_id_pop <-
  cases_raw %>%
  distinct(county, population) %>%
  mutate(population = as.integer(population)) %>%
  bind_rows(
    .,
    left_join(., county_region_key, by = "county") %>%
      select(-county) %>%
      group_by(region) %>%
      summarize(
        population = sum(population),
        .groups = "drop"
      ) %>%
      rename(county = region),
    select(., -county) %>%
      summarize(population = sum(population)) %>%
      mutate(county = "California")
  ) %>%
  mutate(id = seq_len(n()), .before = 1)


cases <-
  cases_raw %>%
  select(county, date, cases, deaths, cumulative_deaths) %>%
  drop_na()

# Hospitalized should be hospitalized_covid_confirmed_patients
# ICU should be icu_covid_confirmed_patients
hosp <-
  read_csv(hosp_url) %>%
  select(county,
    date = todays_date,
    hospitalized = hospitalized_covid_confirmed_patients,
    icu = icu_covid_confirmed_patients
  )

full_dat <-
  full_join(cases, hosp) %>%
  replace_na(list(
    hospitalized = 0,
    icu = 0
  ))

earliest_date_elligible <- ymd("2022-06-01")
latest_date <- max(cases$date, na.rm = TRUE)
last_date_to_report <- latest_date - 2
first_date_to_report <- earliest_date_elligible + (as.numeric(last_date_to_report - earliest_date_elligible) %% time_interval_in_days) + 1

# Should look into re-estimating case and death delay distributions

dat <-
  full_dat %>%
  filter(
    date >= first_date_to_report,
    date <= last_date_to_report
  ) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(
    cases_est_prop_reported = case_reporting_delay_ecdf(days_ago),
    death_est_prop_reported = death_reporting_delay_ecdf(days_ago)
  ) %>%
  mutate(
    est_cases = cases / cases_est_prop_reported,
    est_deaths = deaths / death_est_prop_reported
  ) %>%
  mutate(time = floor(as.numeric(date - min(date)) / time_interval_in_days)) %>%
  group_by(county, time) %>%
  summarize(
    start_date = min(date),
    end_date = max(date),
    cases = sum(cases),
    deaths = sum(deaths),
    est_cases = round(sum(est_cases)),
    est_deaths = round(sum(est_deaths)),
    cumulative_deaths = last(cumulative_deaths),
    hospitalized = last(hospitalized),
    icu = last(icu),
    .groups = "drop"
  ) %>%
  bind_rows(
    .,
    left_join(., county_region_key, by = "county") %>%
      select(-county) %>%
      group_by(region, start_date, end_date, time) %>%
      summarize(across(everything(), sum), .groups = "drop") %>%
      rename(county = region),
    select(., -county) %>%
      group_by(start_date, end_date, time) %>%
      summarize(across(everything(), sum), .groups = "drop") %>%
      mutate(county = "California")
  )

# Initialization Values ---------------------------------------------------

initialization_values <-
  dat %>%
  filter(time == 0) %>%
  select(county,
    H = hospitalized,
    ICU = icu,
    D = cumulative_deaths
  ) %>%
  right_join(county_id_pop, .) %>%
  mutate(across(c(H, ICU, D), ~ . + 1)) %>%
  mutate(remaining_population = population - (H + ICU + D))

write_csv(filter(dat, time != 0), "data/cases_hospitalizations_by_county.csv")
write_csv(initialization_values, "data/initialization_values.csv")
write_csv(county_id_pop, "data/county_id_pop.csv")

# Clear results for next fit
if (dir_exists(results_dir)) {
  dir_delete(results_dir)
}
