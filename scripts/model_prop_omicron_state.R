library(tidyverse)
library(broom)
library(broom.mixed)
library(brms)
library(tidybayes)
library(lubridate)
library(scales)
source("src/plot_functions.R")
options(brms.backend = "cmdstanr",
        mc.cores=parallel::detectCores())

county_variant_data <- read_csv("data/county_greek_data.csv") %>%
  pivot_longer(-c(county, date), names_to = "greek", values_to = "n") %>%
  mutate(greek = fct_other(greek, keep = "omicron", other_level = "other")) %>%
  count(county, date, greek, wt = n) %>%
  pivot_wider(names_from = greek, values_from = n)

state_variant_dat <-
  county_variant_data %>%
  group_by(date) %>%
  summarize(omicron = sum(omicron),
            other = sum(other))


date_scale_info <- scale(state_variant_dat$date)

date_scale_center <- attr(date_scale_info, "scaled:center")
date_scale_scale <- attr(date_scale_info, "scaled:scale")


state_variant_data_for_brms <-
  state_variant_dat %>%
  mutate(date_scaled = (as.numeric(date) - date_scale_center) / date_scale_scale) %>%
  mutate(n = omicron + other)


placeholder_data <-
  crossing(date = seq(min(state_variant_dat$date), today(), 1),
           omicron = 0,
           n = 10) %>%
  mutate(date_scaled = (as.numeric(date) - date_scale_center) / date_scale_scale)


brm_res_1 <- brm(formula = omicron | trials(n) ~ 1 + date_scaled,
                 data = state_variant_data_for_brms,
                 family = binomial)

augment_with_epred_data <-function(model) {
  epred_draws(model, placeholder_data) %>%
    mutate(.epred = .epred / n) %>%
    ungroup() %>%
    select(date, .epred) %>%
    group_by(date) %>%
    summarize(epred = mean(.epred), .groups = "drop")
}

all_epred_data <- augment_with_epred_data(brm_res_1)

prop_omicron_estimates_state <-
  ggplot(mapping = aes(date, epred)) +
  geom_point(data = state_variant_data_for_brms %>%
               mutate(epred = omicron / n),
             mapping = aes(size = n), alpha = 0.5) +
  geom_line(data = all_epred_data) +
  scale_y_continuous(name = "Proportion Omicron", labels = percent) +
  scale_size_continuous(name = "# Sequences", labels = comma) +
  scale_x_date(name = "Date") +
  cowplot::theme_cowplot() +
  theme(legend.position = "bottom") +
  ggtitle("Introduction of Omicron in CA")


save_plot_target_asp(filename = "figures/prop_omicron_estimates_state.pdf", plot = prop_omicron_estimates_state, base_asp = 16/9, base_height = 5)
