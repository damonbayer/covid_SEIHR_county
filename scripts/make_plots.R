library(tidyverse)
library(gridExtra)
library(tidybayes)
source("src/plot_functions.R")

raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases)

county_region_key <- read_csv("data/county_region_key.csv")

generated_quantities_summary <- read_csv("results/generated_quantities_summary.csv")

prior_generated_quantities_summary <-
  read_csv("results/prior_generated_quantities_summary.csv") %>%
  filter(is.na(date)) %>%
  mutate(distribution = "Prior")

all_pp_summaries <- bind_rows(
  read_csv("results/posterior_predictive_summary_CA.csv") %>%
    rename(place_name = county) %>%
    mutate(place_name = "California",
           place_type = "State"),
  read_csv("results/posterior_predictive_summary_regions.csv") %>%
    rename(place_name = region) %>%
    mutate(place_type = "Region"),
  read_csv("results/posterior_predictive_summary_counties.csv") %>%
    rename(place_name = county) %>%
    mutate(place_type = "County"))

dat_tidy <-
  raw_dat %>%
  filter(time > 0) %>%
  select(county, date, est_other_cases, est_omicron_cases, hospitalizations, est_deaths, icu) %>%
  rename(est_death = est_deaths) %>%
  pivot_longer(-c(date, county)) %>%
    left_join(county_region_key) %>%
    rename(place_name = county) %>%
    mutate(place_type = "county") %>%
  bind_rows(.,
            group_by(., date, name) %>%
              summarize(value = sum(value),
                        .groups = "drop") %>%
              mutate(place_name = "California",
                     place_type = "state"),
            group_by(., date, region, name) %>%
              summarize(value = sum(value),
                        .groups = "drop") %>%
              rename(place_name = region) %>%
              mutate(place_type = "region")
            ) %>%
    select(-region)

make_post_pred_plot <- function(target_place_name) {
  tmp_posterior_predictive_intervals <-
    all_pp_summaries %>%
    filter(place_name == target_place_name)

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(place_name == target_place_name)

  target_place_type <- tmp_dat_tidy$place_type[1]

  ggplot(mapping = aes(date, value)) +
    facet_wrap(. ~ name,
               scales = "free_y",
               labeller = as_labeller(
                 c(hospitalizations = "Concurrent Hospitalizations",
                   est_omicron_cases = "Reported Omicron Cases (7 day bins)",
                   est_other_cases = "Reported Other Cases (7 day bins)",
                   icu = "Concurrent ICU",
                   est_death = "Reported Deaths"))) +
    geom_lineribbon(data = tmp_posterior_predictive_intervals,
                    mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_dat_tidy) +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place_name, target_place_type, sep = " ") %>% str_to_title(),
            subtitle = str_c("Forecasted ", max(dat_tidy$date))) +
    my_theme
}

make_time_varying_gq_plot <- function(target_place_name) {
  generated_quantities_summary %>%
    filter(county == target_place_name) %>%
    filter(!is.na(date)) %>%
    ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
    facet_wrap(. ~ name, scales = "free_y") +
    geom_lineribbon() +
    scale_y_continuous(label = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place_name %>% str_to_title(), "County", sep = " "),
            subtitle = str_c("Forecasted", max(dat_tidy$date), sep = " ")) +
    my_theme
}

make_scalar_gq_plot <- function(target_place_name) {
  generated_quantities_summary %>%
    filter(county == target_place_name) %>%
    mutate(distribution = "Posterior") %>%
    bind_rows(prior_generated_quantities_summary) %>%
    filter(is.na(date)) %>%
    ggplot(aes(x = value, y = distribution, xmin = .lower, xmax = .upper, color = distribution)) +
    facet_wrap(. ~ name, scales = "free_x") +
    geom_pointinterval() +
    scale_y_discrete(name = "Distribution") +
    scale_x_continuous(name = "Value") +
    scale_color_discrete(name = "Distribution") +
    theme(legend.position = "bottom")
}


ggsave2(filename = "figures/pp_plots.pdf",
        plot = dat_tidy %>%
          mutate(place_type = fct_relevel(place_type, "state", "region", "county")) %>%
          arrange(place_type, place_name) %>%
          pull(place_name) %>%
          unique() %>%
          map(make_post_pred_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)

ggsave2(filename = "figures/time_varying_gq_plots.pdf",
        plot = dat_tidy %>%
          filter(place_type == "county") %>%
          distinct(place_name) %>%
          arrange(place_name) %>%
          pull(place_name) %>%
          map(make_time_varying_gq_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)

ggsave2(filename = "figures/scalar_gq_plots.pdf",
        plot = dat_tidy %>%
          filter(place_type == "county") %>%
          distinct(place_name) %>%
          arrange(place_name) %>%
          pull(place_name) %>%
          map(make_scalar_gq_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)


# Save only Target Counties -----------------------------------------------
target_counties <- c("Los Angeles", "San Francisco", "Sacramento")

hosp_forecast_target_counties <-
  ggplot(mapping = aes(date, value)) +
  facet_wrap(. ~ place_name,
             scales = "free_y") +
    geom_lineribbon(data = all_pp_summaries %>%
                      filter(place_name %in% target_counties) %>%
                      filter(name == "hospitalizations"),
                    mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = dat_tidy %>%
                 filter(place_name %in% target_counties) %>%
                 filter(name == "hospitalizations")) +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c("County", "Hospitalizations", sep = " ") %>% str_to_title(),
            subtitle = str_c("Forecasted ", max(dat_tidy$date))) +
    my_theme

hosp_forecast_target_regions <-
  ggplot(mapping = aes(date, value)) +
  facet_wrap(. ~ place_name,
             scales = "free_y") +
  geom_lineribbon(data = all_pp_summaries %>%
                    filter(place_type == "Region") %>%
                    filter(name == "hospitalizations"),
                  mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = dat_tidy %>%
               filter(place_type == "region") %>%
               filter(name == "hospitalizations")) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_date(name = "Date") +
  ggtitle(str_c("Regional", "Hospitalizations", sep = " ") %>% str_to_title(),
          subtitle = str_c("Forecasted ", max(dat_tidy$date))) +
  my_theme


hosp_forecast_target_state <-
  ggplot(mapping = aes(date, value)) +

  geom_lineribbon(data = all_pp_summaries %>%
                    filter(place_type == "State") %>%
                    filter(name == "hospitalizations"),
                  mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = dat_tidy %>%
               filter(place_type == "state") %>%
               filter(name == "hospitalizations")) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_date(name = "Date") +
  ggtitle(str_c("State", "Hospitalizations", sep = " ") %>% str_to_title(),
          subtitle = str_c("Forecasted ", max(dat_tidy$date))) +
  my_theme


save_plot_target_asp(filename = "figures/hosp_forecast_target_counties.pdf",
                     plot = hosp_forecast_target_counties,
                     ncol = 3,
                     base_asp = 16/9, base_height = 6)

save_plot_target_asp(filename = "figures/hosp_forecast_target_regions.pdf",
                     plot = hosp_forecast_target_regions,
                     ncol = 3, nrow = 2,
                     base_asp = 16/9, base_height = 3)

save_plot_target_asp(filename = "figures/hosp_forecast_target_state.pdf",
                     plot = hosp_forecast_target_state,
                     base_asp = 16/9, base_height = 4)
