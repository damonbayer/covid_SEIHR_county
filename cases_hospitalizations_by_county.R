library(tidyverse)
library(ckanr)
library(lubridate)

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

# case_reporting_delay_ecdf <- read_rds("/Users/damon/Documents/covid_SEIHR_county/case_reporting_delay_ecdf.rds")
#
# case_reporting_delay_ecdf(7)
#
#
# tibble(delay = 1:30,
#        prop_reported = case_reporting_delay_ecdf(delay)) %>%
#   ggplot(aes(delay, prop_reported)) +
#   geom_line() +
#   geom_point() +
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(breaks = seq(0,28, by = 7)) +
#   cowplot::theme_minimal_grid()


full_dat <- full_join(cases, hosp) %>%
  select(date, county, cases, hospitalized_covid_patients)

dat <-
  full_dat %>%
  drop_na() %>%
  filter(date >= "2021-11-28") %>%
  mutate(days_ago = as.numeric(max(date) - date)) %>%
  filter(days_ago >= 6) %>%
  filter(date >= days_ago[which(as.numeric(date - max(date)) %% 7 == 0)[1]]) %>%
  mutate(time = floor(as.numeric(date - min(date)) / 7)) %>%
  group_by(time, county) %>%
  summarize(date = last(date),
            cases = sum(cases),
            hospitalizations = last(hospitalized_covid_patients), .groups = "drop")


county_id_key <-
  dat %>%
  select(county) %>%
  distinct() %>%
  mutate(id = 1:n()) %>%
  select(id, county)

write_csv(dat, "data/cases_hospitalizations_by_county.csv")
write_csv(county_id_key, "data/county_id_key.csv")
