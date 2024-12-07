library(tidyverse)


plotData <- function(dat, ydata, ylabel, ...) {
  ggplot(dat, aes(x = indVar, y = {{ydata}}, colour = indVar, fill = indVar)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, na.rm = TRUE) +
    geom_point(shape = 1, position = position_jitter(0.2, seed = 43), na.rm = TRUE) +
    scale_y_continuous(name = ylabel, ...) +
    labs(x = "Individual variation", y = ylabel) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(w ~ model, labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}


dat <- read_rds("../data/data-processed.rds") |>
  mutate(w = paste0("omega == ", w)) |>
  mutate(indVar = fct_relevel(indVar, "low", "medium", "high"))

dat |> plotData(richness, "Species richness")
ggsave("../figures/richness.pdf", width = 6, height = 5.5)

dat |> plotData(diversity, "Species diversity", limits = c(0, NA))
ggsave("../figures/diversity.pdf", width = 6, height = 5.5)

dat |> plotData(clustering, "Trait clustering (square root scale)", trans = "sqrt")
ggsave("../figures/clustering.pdf", width = 6, height = 5.5)

dat |> plotData(robustness, "Average community robustness (log scale)", trans = "log10")
ggsave("../figures/robustness.pdf", width = 6, height = 5.5)
