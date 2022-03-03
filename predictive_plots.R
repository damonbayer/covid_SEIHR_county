library(tidyverse)
library(tidybayes)
prior_predictive_intervals <-
  read_csv("results/prior_predictive.csv") %>%
  pivot_longer(-c(iteration, chain)) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  select(time, name, value) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = time_interval_in_days, length.out = max(.$time)))) %>%
  mutate(source = "Prior")

map_predictive_samples_intervals <-
  read_csv("results/map_predictive.csv") %>%
  pivot_longer(-c(iteration, chain)) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  select(time, name, value) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = time_interval_in_days, length.out = max(.$time)))) %>%
  mutate(source = "MAP")


mle_predictive_samples_intervals <-
  read_csv("results/mle_predictive.csv") %>%
  pivot_longer(-c(iteration, chain)) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  select(time, name, value) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  left_join(.,tibble(time = 1:max(.$time),
                     date = seq(min(raw_dat$date), by = time_interval_in_days, length.out = max(.$time)))) %>%
  mutate(source = "MLE")



map_predictive_samples_intervals %>%
  select(date, time) %>%
  distinct()
max(map_predictive_samples_intervals$time)


ggplot(mapping = aes(date, value)) +
  facet_grid(name ~ source, scale = "free_y") +
  geom_lineribbon(data = bind_rows(prior_predictive_intervals,
                                   map_predictive_samples_intervals,
                                   mle_predictive_samples_intervals) %>%
                    filter(date <= max(dat_tidy$date)),
                  mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = dat_tidy %>% filter(county == "Orange") %>% crossing(source = c("Prior", "MAP", "MLE"))) +
  scale_y_continuous(labels = comma) +
  scale_fill_brewer() +
  cowplot::theme_minimal_grid()
