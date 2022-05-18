library(tidyverse)
library(gridExtra)

raw_dat <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
  rename(new_cases = cases)

county_region_key <- read_csv("data/county_region_key.csv")

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
  select(county, date, est_other_cases, est_omicron_cases, hospitalizations) %>%
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
                   est_omicron_cases = "Reported Omicron Cases (3 day bins)",
                   est_other_cases = "Reported Other Cases (3 day bins)"))) +
    geom_lineribbon(data = tmp_posterior_predictive_intervals,
                    mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_dat_tidy) +
    scale_fill_brewer(name = "Credible Interval Width") +
    scale_y_continuous(name = "Count", labels = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place_name, target_place_type, sep = " ") %>% str_to_title(),
            subtitle = str_c("Forecasted ", max(dat_tidy$date))) +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "bottom")
}


ggsave2(filename = "figures/pp_plots.pdf",
        plot = dat_tidy %>%
          pull(place_name) %>%
          unique() %>%
          map(make_post_pred_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)
