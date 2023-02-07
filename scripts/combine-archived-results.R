library(dplyr)
library(stringr)

# Specify column names for each format, only used for reference
# A new format may have name mismatches i.e. cum_death & cumulative_deaths

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
  "creation_date" = str_split_fixed(name_details, pattern = "-format", n = 2)[, 1],
  "format" = str_split_fixed(name_details, pattern = "-format", n = 2)[, 2]
)

all_forecasts <- readr::read_csv(
  here::here("results-archive", paste0(archived_file_names_df$file_names[1])),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(creation_date = archived_file_names_df$creation_date[1])

for (i in 2:nrow(archived_file_names_df)) {
  single_forecast <- readr::read_csv(
    here::here("results-archive", paste0(archived_file_names_df$file_names[i])),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(creation_date = archived_file_names_df$creation_date[i])

  if (archived_file_names_df$format[i] %in% 3:5) {
    single_forecast <- rename(single_forecast, cumulative_deaths = cum_death)
  }

  all_forecasts <- all_forecasts %>%
    bind_rows(single_forecast)

}

# If you want to save the combined forecasts
readr::write_csv(all_forecasts, here::here("results-archive", "combined_archived_forecasts.csv"))
