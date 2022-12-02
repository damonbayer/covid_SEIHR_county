library(tidyverse)
tmp <- read_csv("notes/county_greek_data.csv")

tmp %>%
  select(county, date, delta, omicron) %>%
  pivot_longer(-c(county, date)) %>%
  count(date, name, wt = value) %>%
  group_by(date) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(date, prop, color = name)) +
  geom_smooth() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = lubridate::ymd("2022-03-01"))

