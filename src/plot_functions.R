library(ggplot2)
library(scales)
library(cowplot)
theme_set(cowplot::theme_minimal_grid())
save_plot_target_asp <- function (filename, plot, ncol = 1, nrow = 1, base_height = 3.71,
                                  base_asp = 1730/650, base_width = NULL) {
  cowplot::save_plot(filename, plot, ncol = ncol, nrow = nrow, base_height = base_height,
                     base_asp = base_asp * nrow / ncol, base_width = base_width)
}


my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))
