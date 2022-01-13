library(gridExtra)
library(tidyverse)
library(scales)
library(tidybayes)
library(fs)
library(cowplot)

raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases)

county_id_key <- read_csv("data/county_id_key.csv")

raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases)

dat_tidy <-
  raw_dat %>%
  filter(time > 0) %>%
  select(county, date, new_cases, hospitalizations) %>%
  pivot_longer(-c(date, county))

posterior_predictive_intervals <-
  dir_ls("//dfs6/pub/bayerd/covid_SEIHR_county/results") %>%
  enframe(name = NULL, value = "full_path") %>%
  mutate(file_name = full_path %>%
           path_file() %>%
           path_ext_remove()) %>%
  filter(str_ends(file_name, "\\d+")) %>%
  mutate(county_id = file_name %>%
           str_extract("\\d+$") %>%
           as.numeric()) %>%
  filter(str_detect(file_name, "posterior_predictive")) %>%
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
  select(county_id, time, name, value) %>%
  group_by(county_id, time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = 7, length.out = max(.$time)))) %>%
  left_join(county_id_key %>% rename(county_id = id)) %>%
  select(county, date, name, value, starts_with("."))

make_plot <- function(county_name) {
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
                   new_cases = "Weekly Cases"))) +
    geom_lineribbon(data = tmp_posterior_predictive_intervals,
                    mapping = aes(date, value, ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_dat_tidy,
               mapping = aes(date, value)) +
    scale_fill_brewer(name = "Credible Interval Width") +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(label = str_c(county_name, " County"),
            subtitle = str_c("Forecasted ", max(tmp_dat_tidy$date) + 7)) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "bottom")
}


plot_tibble <-
  tibble(county_name = unique(posterior_predictive_intervals$county)) %>%
  mutate(plot_obj = map(county_name, make_plot))

ggsave2(filename = "plots.pdf",
        plot = marrangeGrob(plot_tibble$plot_obj, nrow=1, ncol=1),
        width = 11,
        height = 8.5)

prior_gq_samples <- read_csv("results/prior_gq_samples.csv")

prior_plot <-
  prior_gq_samples %>%
  pivot_longer(-c(iteration, chain)) %>%
  ggplot(aes(value)) +
  facet_wrap(. ~ name, scales = "free") +
  stat_halfeye(normalize = "panels") +
  cowplot::theme_cowplot() +
  labs(title = "Model Priors",
       x = NULL,
       y = NULL)

ggsave2(filename = "prior_plot.pdf", plot = prior_plot)