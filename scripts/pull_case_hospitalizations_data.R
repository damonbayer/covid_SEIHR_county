library(tidyverse)
library(ckanr)
library(lubridate)
library(splines)
library(fs)

results_dir <- "results"

case_reporting_delay_ecdf <- read_rds("data/case_reporting_delay_ecdf.rds")
death_reporting_delay_ecdf <- read_rds("data/death_delay_ecdf.rds")
county_region_key <- read_csv("data/county_region_key.csv")

prop_omicron_county_dat <-
  read_csv("data/prop_omicron_county_dat.csv") %>%
  # Repeat last prop if data isn't updated
  bind_rows(., crossing(group_by(., county) %>%
                          filter(date == max(date)) %>%
                          select(-date),
                        date = seq(max(.$date), today(), 1)))

# variants_dat <- read_tsv("https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/master/data/omicron-us/omicron-us_location-variant-sequence-counts.tsv") %>%
#   filter(location == "California") %>%
#   select(-location) %>%
#   distinct() %>%
#   pivot_wider(names_from = variant, values_from = sequences, values_fill = 0) %>%
#   pivot_longer(-date, names_to = "variant", values_to = "sequences") %>%
#   mutate(sequences = sequences + 1) %>%
#   group_by(date) %>%
#   summarize(variant = variant,
#             prop = sequences / sum(sequences)) %>%
#   filter(variant == "Omicron") %>%
#   select(date, prop_omicron = prop)

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


ckanr_setup(url="https://data.ca.gov")
ckan <- quiet(ckanr::src_ckan("https://data.ca.gov"))

# get resources
resources <- rbind(resource_search("name:covid-19", as = "table")$results,
                   resource_search("name:hospitals by county", as = "table")$results)


cases_deaths_url <- resources %>% filter(name == "Statewide COVID-19 Cases Deaths Tests") %>% pull(url)
hosp_url <- resources %>% filter(name == "Statewide Covid-19 Hospital County Data") %>% pull(url)

cases <-
  read_csv(cases_deaths_url) %>%
  filter(area_type == "County") %>%
  mutate(date = lubridate::ymd(date),
         deaths = as.integer(deaths),
         reported_cases = as.integer(reported_cases),
         cases = as.integer(cases),
         positive_tests = as.integer(positive_tests),
         total_tests = as.integer(total_tests)) %>%
  select(date,
         cases = cases,
         tests = total_tests,
         deaths,
         county = area) %>%
  arrange(date, county)

hosp <-
  read_csv(hosp_url) %>%
  mutate(todays_date = lubridate::ymd(todays_date),
         hospitalized_covid_patients = as.integer(hospitalized_covid_confirmed_patients),
         icu_covid_confirmed_patients = as.integer(icu_covid_confirmed_patients),
         icu_suspected_covid_patients = as.integer(icu_suspected_covid_patients)) %>%
  replace_na(list(hospitalized_covid_patients = 0L,
                  icu_covid_confirmed_patients = 0L,
                  icu_suspected_covid_patients = 0L)) %>%
  mutate(icu_covid_patients = icu_covid_confirmed_patients) %>%
  select(date = todays_date,
         hospitalized_covid_patients,
         icu_covid_patients,
         county) %>%
  mutate(hospitalized_covid_patients = hospitalized_covid_patients - icu_covid_patients) %>%
  mutate(hospitalized_covid_patients = ifelse(hospitalized_covid_patients < 0, 0, hospitalized_covid_patients))

county_pop <- read_csv("data/county_pop.csv") %>% rename_with(str_to_lower)

full_dat <-
  full_join(cases, hosp) %>%
  select(date, county, cases, tests, hospitalized_covid_patients, icu_covid_patients, deaths)

time_interval_in_days <- 7

earliest_date_elligible_to_report <- ymd("2021-12-12")
latest_date <- max(full_dat$date, na.rm = T)
last_date_to_report <- latest_date - 2
first_date_to_report <-  earliest_date_elligible_to_report + (as.numeric(last_date_to_report - earliest_date_elligible_to_report) %% time_interval_in_days) + 1

cum_deaths <- cases %>%
  drop_na() %>%
  group_by(county) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(death_est_prop_reported = death_reporting_delay_ecdf(days_ago)) %>%
  mutate(est_deaths = round(deaths/ death_est_prop_reported)) %>%
  mutate(cum_deaths = cumsum(deaths),
         cum_est_deaths = cumsum(est_deaths)) %>%
  ungroup() %>%
  select(date, county, cum_deaths, cum_est_deaths)

dat <-
  full_dat %>%
  drop_na() %>%
  filter(date >= first_date_to_report,
         date <= last_date_to_report) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(est_prop_reported = case_reporting_delay_ecdf(days_ago),
         death_est_prop_reported = death_reporting_delay_ecdf(days_ago)) %>%
  mutate(est_cases = round(cases / est_prop_reported),
         est_tests = round(tests / est_prop_reported),
         est_deaths = round(deaths/ death_est_prop_reported)) %>%
  left_join(prop_omicron_county_dat) %>%
  mutate(est_omicron_cases = prop_omicron_cases * est_cases,
         est_other_cases = (1- prop_omicron_cases) * est_cases,
         est_omicron_tests = prop_omicron_cases * est_tests,
         est_other_tests = (1- prop_omicron_cases) * est_tests) %>%
  mutate(lump = floor(as.numeric(date - min(date)) / time_interval_in_days) + 1) %>%
  group_by(lump, county) %>%
  summarize(time = lump[1] * time_interval_in_days / 7,
            date = last(date),
            cases = sum(cases),
            est_cases = sum(est_cases),
            est_omicron_cases = round(sum(est_omicron_cases)),
            est_other_cases = round(sum(est_other_cases)),
            tests = sum(tests),
            est_tests = sum(est_tests),
            est_omicron_tests = round(sum(est_omicron_tests)),
            est_other_tests = round(sum(est_other_tests)),
            hospitalizations = last(hospitalized_covid_patients),
            deaths = sum(deaths),
            est_deaths = sum(est_deaths),
            icu = last(icu_covid_patients),
            .groups = "drop") %>%
  select(-lump) %>%
  mutate(est_tests = if_else(est_other_tests == 0, est_tests + 1, est_tests),
         est_other_tests = if_else(est_other_tests == 0,  1, est_other_tests)) %>%
  left_join(cum_deaths, by = c("date", "county")) %>%
  mutate(prev_cum_deaths = cum_deaths - deaths,
         prev_cum_est_deaths = cum_est_deaths - est_deaths) %>%
  bind_rows(.,
            left_join(., county_region_key) %>%
              select(-county) %>%
              group_by(region, date, time) %>%
              summarize(across(everything(), sum), .groups = "drop") %>%
              rename(county = region),
            select(., -county) %>%
              group_by(date, time) %>%
              summarize(across(everything(), sum), .groups = "drop") %>%
              mutate(county = "California"))

initialization_values <-
  full_dat %>%
  drop_na() %>%
  filter(date >= first_date_to_report - 6,
         date <= first_date_to_report) %>%
  left_join(prop_omicron_county_dat) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(est_prop_reported = case_reporting_delay_ecdf(days_ago)) %>%
  mutate(est_cases = round(cases / est_prop_reported)) %>%
  mutate(est_omicron_cases = round(prop_omicron_cases * est_cases),
         est_other_cases = round((1 - prop_omicron_cases) * est_cases)) %>%
  group_by(county) %>%
  summarize(est_cases = sum(est_cases),
            est_omicorn_cases = sum(est_omicron_cases),
            est_other_cases = sum(est_other_cases),
            hospitalizations = last(hospitalized_covid_patients),
            icu = last(icu_covid_patients)) %>%
  pivot_longer(-county) %>%
  mutate(value = if_else(value == 0, 1, value)) %>%
  pivot_wider(county) %>%
  bind_rows(.,
            left_join(., county_region_key) %>%
              select(-county) %>%
              group_by(region) %>%
              summarize(across(everything(), sum), .groups = "drop") %>%
              rename(county = region),
            select(., -county) %>%
              summarize(across(everything(), sum), .groups = "drop") %>%
              mutate(county = "California"))

county_id_key <-
  dat %>%
  select(county) %>%
  distinct() %>%
  mutate(id = 1:n()) %>%
  select(id, county)

write_csv(dat, "data/cases_hospitalizations_by_county.csv")
write_csv(initialization_values, "data/initialization_values.csv")
write_csv(county_id_key, "data/county_id_key.csv")
# Clear results for next fit
if (dir_exists(results_dir)) {
  dir_delete(results_dir)
}
