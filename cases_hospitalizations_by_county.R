library(tidyverse)
library(ckanr)
library(lubridate)

case_reporting_delay_ecdf <- read_rds("case_reporting_delay_ecdf.rds")

variants_dat <- read_tsv("https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/master/data/omicron-us/omicron-us_location-variant-sequence-counts.tsv") %>%
  filter(location == "California") %>%
  select(-location) %>%
  group_by(date) %>%
  summarize(variant = variant,
            prop = sequences / sum(sequences)) %>%
  filter(variant == "Omicron") %>%
  select(date, prop_omicron = prop)

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

# resource_ids <- list(cases = resources$resource_id[resources$name == "COVID-19 Cases"],
#                      tests = resources$resource_id[resources$name == "COVID-19 Testing"],
#                      hosp = resources$resource_id[resources$name == "Hospitals By County"])

resource_ids <- list(cases_deaths = resources$id[resources$name == "Statewide COVID-19 Cases Deaths Tests"],
                     hosp = resources$id[resources$name == "Statewide Covid-19 Hospital County Data"])


# pull resources into data frames (adds extra cols _id and _full_text)
cases <-
  tbl(src = ckan$con, from = resource_ids$cases_deaths) %>%
  as_tibble() %>%
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
  tbl(src = ckan$con, from = resource_ids$hosp) %>%
  as_tibble() %>%
  mutate(todays_date = lubridate::ymd(todays_date),
         hospitalized_covid_patients = as.integer(hospitalized_covid_patients),
         icu_covid_confirmed_patients = as.integer(icu_covid_confirmed_patients),
         icu_suspected_covid_patients = as.integer(icu_suspected_covid_patients)) %>%
  mutate(icu_covid_patients =
           if_else(is.na(icu_covid_confirmed_patients), 0L, icu_covid_confirmed_patients) +
           if_else(is.na(icu_suspected_covid_patients), 0L, icu_suspected_covid_patients)) %>%
  select(date = todays_date,
         hospitalized_covid_patients,
         icu_covid_patients,
         county)


county_pop <- read_csv("data/county_pop.csv") %>% rename_all(str_to_lower)

full_dat <- full_join(cases, hosp) %>%
  select(date, county, cases, hospitalized_covid_patients)


latest_date <- max(full_dat$date, na.rm = T) + 1
prop_omicron_model <- glm(prop_omicron ~ bs(date), data = variants_dat, family = gaussian(link = "logit"))

time_interval_in_days <- 3

dat <-
  full_dat %>%
  drop_na() %>%
  filter(date >= "2021-12-12",
         date <= latest_date - 6) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(est_prop_reported = case_reporting_delay_ecdf(days_ago)) %>%
  mutate(est_cases = round(cases / est_prop_reported)) %>%
  mutate(lump = floor(as.numeric(date - min(date)) / time_interval_in_days) + 1) %>%
  group_by(lump, county) %>%
  summarize(time = lump[1] * time_interval_in_days / 7,
            date = last(date),
            cases = sum(cases),
            est_cases = sum(est_cases),
            hospitalizations = last(hospitalized_covid_patients),
            .groups = "drop") %>%
  mutate(.,
         prop_omicron_cases = predict(prop_omicron_model,
                                      newdata = .,
                                      type = "response"),
         prop_omicron_hospitalizations = predict(prop_omicron_model,
                                                 newdata = mutate(., date = date - 7),
                                                 type = "response")) %>%
  select(-lump)

initialization_values <-
  full_dat %>%
  drop_na() %>%
  filter(date >= min(dat$date) - time_interval_in_days - 6,
         date <= min(dat$date) - time_interval_in_days) %>%
  mutate(days_ago = as.numeric(latest_date - date)) %>%
  mutate(est_prop_reported = case_reporting_delay_ecdf(days_ago)) %>%
  mutate(est_cases = round(cases / est_prop_reported)) %>%
  group_by(county) %>%
  summarize(est_cases = sum(est_cases),
            hospitalizations = last(hospitalized_covid_patients))


county_id_key <-
  dat %>%
  select(county) %>%
  distinct() %>%
  mutate(id = 1:n()) %>%
  select(id, county)

write_csv(dat, "data/cases_hospitalizations_by_county.csv")
write_csv(initialization_values, "data/initialization_values.csv")
write_csv(county_id_key, "data/county_id_key.csv")
