dat %>%
  filter(county == "Orange") %>%
  select(county, time, date, starts_with("est")) %>%
  # mutate(est_pos = est_cases / est_tests,
  #        est_omicron_pos = est_omicron_cases / est_omicron_tests,
  #        est_other_pos = est_other_cases / est_other_tests) %>%
  pivot_longer(-c(county, time, date)) %>%
  mutate(panel = str_extract(name, "(?<=_)[:alpha:]+$"),
         type = str_extract(name, "(?<=est_)[:alpha:]+(?=_|$)")) %>%
  mutate(type = if_else(type == "tests" | type == "cases", "total", type)) %>%
  ggplot(aes(date, value, group = name, color = type)) +
  facet_wrap(. ~ panel, scale = "free_y") +
  geom_line() +
  geom_point() +
  cowplot::theme_minimal_grid() +
  scale_y_continuous(labels = comma)


variants_dat %>%
  filter(date >= ymd("2021-12-15")) %>%
  ggplot(aes(date, prop_omicron)) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ bs(x),  method.args = list(family = gaussian(link = "logit")), se = F) +
  cowplot::theme_minimal_grid() +
  scale_y_continuous(labels = percent)
tmp <-
  posterior_gq_samples %>%
  select(-county) %>%
  filter(name %in% c("ρ_omicron", "ρ_non_omi")) %>%
  group_by(name) %>%
  median_qi() %>%
  select(name, value, .lower, .upper)


scales::scientific(x = 22)

tmp %>%
  mutate(across(where(is.numeric), ~scales::scientific(x = .)))
