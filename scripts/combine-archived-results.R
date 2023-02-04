library(stringr)

# specify column names for each format
format1 <- c(
  "date", "county", "quantile", "hosp_census_with_covid", "cases"
)
format2 <- c(
  "date", "county", "quantile", "hosp_census_with_covid", "cases",
  "icu", "death"
)
format3 <- c(
  "date", "county", "quantile", "hosp_census_with_covid", "cases",
  "icu", "death", "prev_cum_est_deaths", "cum_death"
)
format4 <- c(
  "date", "county", "quantile", "hosp_census_with_covid", "cases",
  "icu", "death", "cum_death", "prev_cum_est_deaths", "total_cum_death"
)
format5 <- c(
  "date", "county", "quantile", "hosp_census_with_covid", "cases",
  "icu", "death", "cum_death", "prev_cum_est_deaths", "total_cum_death", "Rt"
)
format6 <- c(
  "county", "date", "quantile", "hosp_census_with_covid", "icu", "cases",
  "deaths", "cumulative_deaths", "Rt"
)

# read in each csv and add date column
archived_file_names <- list.files(
  here::here("results-archive"),
  pattern = "results_calcat_format_"
)

name_details <- fs::path_ext_remove(str_split_fixed(
  archived_file_names,
  pattern = "results_calcat_format_",
  n = 2
)[, 2])

archived_file_names_df <- data.frame(
  "file_names" = archived_file_names,
  "update_date" = str_split_fixed(name_details, pattern = "-format", n = 2)[, 1],
  "format" = str_split_fixed(name_details, pattern = "-format", n = 2)[, 2]
)


single_forecast <- readr::read_csv(
  here::here("results-archive", "results_calcat_format_2022-05-31-format1.csv"),
  show_col_types = FALSE)

all_forecasts <- single_forecast

# combine data frames to have most one data frame for all forecasts

