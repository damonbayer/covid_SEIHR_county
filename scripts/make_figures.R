library(tidyverse)
library(gridExtra)
library(tidybayes)
library(fs)
source("src/plot_functions.R")

raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv")

county_region_key <- read_csv("data/county_region_key.csv")

place_type_key <-
  raw_dat %>%
  distinct(county) %>%
  mutate(place_type =
           case_when(
             county %in% county_region_key$county ~ "county",
             county %in% county_region_key$region ~ "region",
             TRUE ~ "state")) %>%
  deframe()

dat_tidy <-
  raw_dat %>%
  select(county,
         date = end_date,
         new_cases = est_cases,
         new_deaths = est_deaths,
         hospitalizations = hospitalized,
         icu) %>%
  pivot_longer(-c(county,date))

results_dir <- "results"

predictive_summary <-
  bind_rows(
    read_csv(path(results_dir, "tidy_posterior_predictive", ext = "csv")) %>%
      mutate(distribution = "posterior"),
    read_csv(path(results_dir, "tidy_prior_predictive", ext = "csv")) %>%
      mutate(distribution = "prior")
  )

scalar_generated_quantities_summary <-
  bind_rows(
    read_csv(path(results_dir, "tidy_scalar_posterior_generated_quantities", ext = "csv")) %>%
      mutate(distribution = "posterior"),
    read_csv(path(results_dir, "tidy_scalar_prior_generated_quantities", ext = "csv")) %>%
      mutate(distribution = "prior")
  ) %>%
  mutate(name = case_when(
    name == "σ_R₀" ~ "sigma_R0",
    name == "ϕ_cases" ~ "phi_cases",
    name == "ϕ_deaths" ~ "phi_deaths",
    name == "ϕ_hosp" ~ "phi_hosp",
    name == "ϕ_icu" ~ "phi_icu",
    TRUE ~ name))

vector_generated_quantities_summary <-
  bind_rows(
    read_csv(path(results_dir, "tidy_vector_posterior_generated_quantities", ext = "csv")) %>%
      mutate(distribution = "posterior"),
    read_csv(path(results_dir, "tidy_vector_prior_generated_quantities", ext = "csv")) %>%
      mutate(distribution = "prior")
  ) %>%
  mutate(name = case_when(
    name == "R₀_t" ~ "R0_t",
    name == "Rₜ_t" ~ "Rt_t",
    name == "β_t" ~ "beta_t",
    TRUE ~ name))

make_posterior_predictive_plot <- function(target_place) {
  ggplot(mapping = aes(date, value)) +
    facet_wrap(~name, scales = "free_y") +
    geom_lineribbon(data = predictive_summary %>%
                      filter(distribution == "posterior") %>%
                      filter(county == target_place),
                    mapping = aes(ymin = .lower,
                                  ymax = .upper)
    ) +
    geom_point(data = dat_tidy %>%
                 filter(county == target_place)) +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place, place_type_key[target_place], sep = " ") %>% str_to_title(),
            subtitle = str_c("Data Through ", max(dat_tidy$date))) +
    my_theme
}

make_vector_generated_quantities_plot <- function(target_place) {
  vector_generated_quantities_summary %>%
    filter(distribution == "posterior") %>%
    filter(county == target_place) %>%
    filter(!is.na(date)) %>%
    ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
    facet_wrap(. ~ name, scales = "free_y") +
    geom_lineribbon() +
    scale_y_continuous(name = "Value", label = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place, place_type_key[target_place], sep = " ") %>% str_to_title(),
            subtitle = str_c("Data Through ", max(dat_tidy$date))) +
    my_theme
}

make_scalar_generated_quantities_plot <- function(target_place) {
  scalar_generated_quantities_summary %>%
    filter(county == target_place) %>%
    mutate(distribution = str_to_title(distribution)) %>%
    ggplot(aes(x = value, y = distribution, xmin = .lower, xmax = .upper, color = distribution)) +
    facet_wrap(. ~ name, scales = "free_x") +
    geom_pointinterval() +
    scale_y_discrete(name = "Distribution") +
    scale_x_continuous(name = "Value") +
    scale_color_discrete(name = "Distribution") +
    theme(legend.position = "bottom") +
    ggtitle(str_c(target_place, place_type_key[target_place], sep = " ") %>% str_to_title(),
            subtitle = str_c("Data Through ", max(dat_tidy$date)))
}


# Save figures ------------------------------------------------------------
dir_create("figures")

ggsave2(filename = path("figures", "posterior_predictive_plots", ext = "pdf"),
        plot = place_type_key %>%
          names() %>%
          map(make_posterior_predictive_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 16,
        height = 9)

ggsave2(filename = path("figures", "vector_generated_quantities_plots", ext = "pdf"),
        plot = place_type_key %>%
          names() %>%
          map(make_vector_generated_quantities_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 16,
        height = 9)

ggsave2(filename = path("figures", "scalar_generated_quantities_plots", ext = "pdf"),
        plot = place_type_key %>%
          names() %>%
          map(make_scalar_generated_quantities_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 16,
        height = 9)
