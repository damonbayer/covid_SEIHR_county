library(gridExtra)
library(tidyverse)
library(scales)
library(tidybayes)
library(fs)
library(cowplot)
library(lubridate)

if (Sys.info()[["sysname"]] == "Linux") {
  results_dir <- "//dfs6/pub/bayerd/covid_SEIHR_county/results"
} else if (Sys.info()[["sysname"]] == "Darwin") {
  results_dir <- "results"
}

time_interval_in_days <- 7


# Process Data ------------------------------------------------------------
raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases)

county_id_key <- read_csv("data/county_id_key.csv") %>%
  add_row(id = 0, county = "California") %>%
  arrange(id)

dat_tidy <-
  raw_dat %>%
  filter(time > 0) %>%
  select(county, date, est_other_cases, est_omicron_cases, hospitalizations) %>%
  pivot_longer(-c(date, county)) %>%
  bind_rows(., group_by(., date, name) %>%
              summarize(value = sum(value),
                        .groups = "drop") %>%
              mutate(county = "California"))


# Predictive Distributions ------------------------------------------------
posterior_predictive_samples <-
  path(results_dir, "posterior_predictive") %>%
  dir_ls() %>%
  enframe(name = NULL, value = "full_path") %>%
  mutate(file_name = full_path %>%
           path_file() %>%
           path_ext_remove()) %>%
  mutate(county_id = file_name %>%
           str_extract("(?<=county_id=)\\d+") %>%
           as.numeric()) %>%
  mutate(results = full_path %>%
           map(read_csv)) %>%
  select(county_id, results) %>%
  unnest(results) %>%
  pivot_longer(-c(county_id, iteration, chain)) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  bind_rows(., group_by(., chain, iteration, time, name) %>%
              summarize(value = sum(value),
                        .groups = "drop") %>%
              mutate(county_id = 0)) %>%
  arrange(county_id)

posterior_predictive_intervals <-
  posterior_predictive_samples %>%
  select(county_id, time, name, value) %>%
  drop_na() %>%
  group_by(county_id, time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = time_interval_in_days, length.out = max(.$time)))) %>%
  left_join(county_id_key %>% rename(county_id = id)) %>%
  select(county, date, name, value, starts_with("."))

posterior_predictive_LEMMA_format <-
  posterior_predictive_samples %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(cases = est_other_cases + est_omicron_cases) %>%
  select(-starts_with("est")) %>%
  select(county_id, time, hospitalizations, cases) %>%
  pivot_longer(c(hospitalizations, cases)) %>%
  group_by(county_id, time, name) %>%
  median_qi(.width = 0.9) %>%
  pivot_longer(cols = c(value, .lower, .upper), names_to = "value_type") %>%
  mutate(quantile = 0.5 +
           case_when(
             value_type == ".lower" ~ -1,
             value_type == "value" ~ -0,
             value_type == ".upper" ~ 1) *
           .width / 2 ) %>%
  select(county_id, time, quantile, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = time_interval_in_days, length.out = max(.$time)))) %>%
  left_join(county_id_key %>% rename(county_id = id)) %>%
  select(date, county, quantile, hosp_census_with_covid = hospitalizations, cases) %>%
  arrange(county, quantile, date)

rm(posterior_predictive_samples)

write_csv(posterior_predictive_LEMMA_format, "posterior_predictive_LEMMA_format.csv")

# Generated Quantities ----------------------------------------------------
prior_gq_samples_all <-
  read_csv(path(results_dir, "prior_generated_quantities.csv")) %>%
  pivot_longer(-c(iteration, chain)) %>%
  select(name, value) %>%
  mutate(county = "Prior",
         source = "Prior")

prior_gq_samples <-
  prior_gq_samples_all %>%
  filter(str_detect(name, "\\[\\d+\\]", negate = T))

prior_gq_samples_time_varying <-
  prior_gq_samples_all %>%
  filter(str_detect(name, "\\[\\d+\\]")) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  mutate(name = fct_inorder(name)) %>%
  group_by(county, name, source, time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 0:max(.$time),
                     date = seq(min(raw_dat$date) - time_interval_in_days, by = time_interval_in_days, length.out = max(.$time) + 1)))

rm(prior_gq_samples_all)

posterior_gq_samples_all <-
  path(results_dir, "generated_quantities") %>%
  dir_ls() %>%
  enframe(name = NULL, value = "full_path") %>%
  mutate(file_name = full_path %>%
           path_file() %>%
           path_ext_remove()) %>%
  mutate(county_id = file_name %>%
           str_extract("(?<=county_id=)\\d+") %>%
           as.numeric()) %>%
  mutate(results = full_path %>%
           map(read_csv)) %>%
  select(county_id, results) %>%
  unnest(results) %>%
  pivot_longer(-c(county_id, iteration, chain)) %>%
  left_join(county_id_key %>% rename(county_id = id)) %>%
  select(county, name, value)
#maybe put processing inside map to save ram?
posterior_gq_samples <-
  posterior_gq_samples_all %>%
  filter(str_detect(name, "\\[\\d+\\]", negate = T))

posterior_gq_samples_time_varying <-
  posterior_gq_samples_all %>%
  filter(str_detect(name, "\\[\\d+\\]")) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  group_by(county, name, time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 0:max(.$time),
                     date = seq(min(raw_dat$date) - time_interval_in_days, by = time_interval_in_days, length.out = max(.$time) + 1)))

rm(posterior_gq_samples_all)

# Plot Functions ----------------------------------------------------------

# prior_gq_samples_time_varying %>%
#   ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
#   facet_wrap(. ~ name, scales = "free_y") +
#   geom_lineribbon() +
#   scale_y_continuous(labels = comma) +
#   scale_fill_brewer() +
#   cowplot::theme_minimal_grid()


make_post_pred_plot <- function(county_name) {
  tmp_posterior_predictive_intervals <-
    posterior_predictive_intervals %>%
    filter(county == county_name) %>%
    mutate(.width = percent(.width))

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(county == county_name)

  ggplot() +
    facet_wrap(. ~ name,
               scales = "free_y",
               labeller = as_labeller(
                 c(hospitalizations = "Concurrent Hospitalizations",
                   est_omicron_cases = "Reported Omicron Cases (3 day bins)",
                   est_other_cases = "Reported Other Cases (3 day bins)"))) +
    geom_lineribbon(data = tmp_posterior_predictive_intervals,
                    mapping = aes(date, value, ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_dat_tidy,
               mapping = aes(date, value)) +
    scale_fill_brewer(name = "Credible Interval Width") +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(label = if_else(county_name == "California",
                            county_name,
                            str_c(county_name, " County")),
            subtitle = str_c("Forecasted ", max(dat_tidy$date) + 6)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "bottom")
}

make_prior_post_plot <- function(county_name) {
  posterior_gq_samples %>%
    filter(county == county_name) %>%
    mutate(source = "Posterior") %>%
    bind_rows(prior_gq_samples) %>%
    ggplot(aes(value, group = source, fill = source, color = source)) +
    facet_wrap(. ~ name, scales = "free") +
    stat_halfeye(normalize = "panels", alpha = 0.5) +
    theme_cowplot() +
    labs(x = NULL,
         y = NULL,
         color = "Source",
         fill = "Source",
         title = if_else(county_name == "California",
                         county_name,
                         str_c(county_name, " County")),
         subtitle = str_c("Fit ", today())) +
    theme(legend.position = "bottom")
}

make_time_varying_plot <- function(county_name) {
  if (county_name == "California") {
    return(ggplot())
  }
  posterior_gq_samples_time_varying %>%
    filter(county == county_name) %>%
    ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
    facet_wrap(. ~ name,
               scales = "free_y") +
    geom_lineribbon() +
    scale_fill_brewer(name = "Credible Interval Width") +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(label = if_else(county_name == "California",
                            county_name,
                            str_c(county_name, " County")),
            subtitle = str_c("Forecasted ", max(dat_tidy$date) + 6)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "bottom")
}

# Create and Save Plots ---------------------------------------------------
plot_tibble <-
  tibble(county_name = unique(posterior_predictive_intervals$county)) %>%
  mutate(post_pred_plot_obj = map(county_name, make_post_pred_plot),
         prior_post_plot_obj = map(county_name, make_prior_post_plot),
         time_varying_plot_obj = map(county_name, make_time_varying_plot))

ggsave2(filename = "post_pred_plots.pdf",
        plot = marrangeGrob(plot_tibble$post_pred_plot_obj, nrow=1, ncol=1),
        width = 18,
        height = 6)

ggsave2(filename = "prior_post_plots.pdf",
        plot = marrangeGrob(plot_tibble$prior_post_plot_obj, nrow=1, ncol=1),
        width = 12,
        height = 8)

ggsave2(filename = "time_varying_plots.pdf",
        plot = marrangeGrob(plot_tibble$time_varying_plot_obj[-1], nrow=1, ncol=1),
        width = 12,
        height = 8)
