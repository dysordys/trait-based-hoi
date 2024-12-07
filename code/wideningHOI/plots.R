library(tidyverse)


read_rds("../../data/widen-processed.rds") |>
  summarize(mean = mean(diversity),
            lower = min(diversity),
            upper = max(diversity),
            .by = c(W, model)) |>
  ggplot(aes(x = W, y = mean, ymin = lower, ymax = upper)) +
  geom_line(color = "steelblue") +
  geom_ribbon(fill = "steelblue", alpha = 0.2) +
  labs(x = "Extra width to HOI kernel", y = "Species diversity") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~ model, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("../../figures/widening.pdf", width = 6, height = 3)
