library(tidyverse)
library(tidybayes)
library(fs)
source("src/plot_functions.R")
prev_model_ci <- read_csv("notes/prev_model_ci.csv")

prev_model_ci %>%
  select(county, date, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(D_total = D_non_omicron + D_omicron,
         E_total = E_non_omicron + E_omicron,
         H_total = H_non_omicron + H_omicron,
         I_total = I_non_omicron + I_omicron,
         ICU_total = ICU_non_omicron + ICU_omicron,
         S_total = S_both + S_omicron_only) %>%
  pivot_longer(-c(county, date)) %>%
  filter(!(name %in%
           c("D_non_omicron",
             "D_omicron",
             "D_total",
             "H_non_omicron",
             "H_omicron",
             "H_total",
             "ICU_non_omicron",
             "ICU_omicron",
             "ICU_total"))) %>%
  filter(date >= "2022-06-01") %>%
  ggplot(aes(date, value, group = county)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line() +
  scale_y_continuous(labels = percent)

# Conclusion:
# Start on June 1
# S: 80% - 95%
# E: 0.1% - 0.4%
# I: 0.25% - 2%
# R: 5% - 20%
# D: Data
# H: Data
# ICU: Data
