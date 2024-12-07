library(tidyverse)

smoothstep <- function(n) ifelse(n < 1, (n > 0) * (n^3*(10+n*(-15+6*n))), 1)

tibble(x = seq(-0.25, 1.25, l=501)) %>%
  mutate(q = smoothstep(x)) %>%
  ggplot(aes(x = x, y = q)) +
  geom_line(colour = "steelblue") +
  scale_x_continuous(name = expression(paste(x))) +
  scale_y_continuous(name = expression(paste(Q(x))), breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("../../figures/smoothstep.pdf", width = 4, height = 3)
