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


date_scale_info <- scale(county_variant_data$date)

date_scale_center <- attr(date_scale_info, "scaled:center")
date_scale_scale <- attr(date_scale_info, "scaled:scale")

county_variant_data_for_brms <-
  county_variant_data %>%
  mutate(date_scaled = (as.numeric(date) - date_scale_center) / date_scale_scale) %>%
  mutate(n = omicron + other)

placeholder_data <-
  crossing(county = unique(county_variant_data$county),
           date = seq(min(county_variant_data$date), today(), 1),
           omicron = 0,
           n = 10) %>%
  mutate(date_scaled = (as.numeric(date) - date_scale_center) / date_scale_scale)

brm_res_1 <- brm(formula = omicron | trials(n) ~ (1 | county) + date_scaled,
                     data = county_variant_data_for_brms,
                     family = binomial)

brm_res_2 <- brm(formula = omicron | trials(n) ~ date_scaled | county,
                 data = county_variant_data_for_brms,
                 family = binomial)


augment_with_epred_data <- function(model) {
  epred_draws(model, placeholder_data) %>%
    mutate(.epred = .epred / n) %>%
    ungroup() %>%
    select(date, county, .epred) %>%
    group_by(date, county) %>%
    summarize(epred = mean(.epred), .groups = "drop")
}

all_epred_data <-
  bind_rows(
    augment_with_epred_data(brm_res_1) %>%
      mutate(model = "random intercept"),
    augment_with_epred_data(brm_res_2) %>%
      mutate(model = "random slope and intercept")) %>%
  unite(col = unique_id, county, model, remove = F)


ggplot(all_epred_data, aes(date, epred, group = unique_id))  +
  facet_wrap(. ~ model) +
  geom_line() +
  scale_y_continuous(name = "Proportion Omicron", labels = percent) +
  scale_x_date(name = "Date") +
  cowplot::theme_cowplot()

prop_omicron_estiamtes_county <-
  ggplot(mapping = aes(date, epred)) +
  facet_wrap(. ~ county) +
  geom_point(data = county_variant_data_for_brms %>%
               mutate(epred = omicron / n),
             alpha = 0.5) +
  geom_line(data = all_epred_data, mapping = aes(group = unique_id, color = model)) +
  scale_y_continuous(name = "Proportion Omicron", labels = percent) +
  scale_x_date(name = "Date") +
  cowplot::theme_cowplot() +
  theme(legend.position = "bottom")


save_plot_target_asp(filename = "figures/prop_omicron_estiamtes_county.pdf", plot = prop_omicron_estiamtes_county.pdf, ncol = 8, nrow = 8, base_height = 1.2, base_asp = 16/9)

prop_omicron_county_dat <-
  all_epred_data %>%
  filter(model == "random intercept") %>%
  select(date, county, prop_omicron_cases = epred) %>%
  bind_rows(.,
            filter(., county %in% c("Siskiyou", "Lassen", "Shasta")) %>%
            group_by(date) %>%
            summarize(prop_omicron_cases = mean(prop_omicron_cases)) %>%
            mutate(county = "Modoc"))


write_csv(prop_omicron_county_dat, "data/prop_omicron_county_dat.csv")
