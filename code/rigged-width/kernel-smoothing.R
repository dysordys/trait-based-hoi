library(tidyverse)

hierKernel <- \(x) pnorm(sqrt(2)*x)
w <- 0.5
W <- 0.01
kappa = 1

tab1 <- tibble(z = seq(-1.1, 0.8, l = 501)) |>
  mutate(Pairwise = hierKernel(z/w),
         `Higher-order` = kappa*hierKernel(z/W)) |>
  mutate(Total = Pairwise + `Higher-order`) |>
  pivot_longer(!z) |>
  mutate(name = fct_relevel(name, "Total", "Pairwise"),
         type = "No trait variance")

tab2 <- tibble(z = seq(-1.1, 0.8, l = 501)) |>
  mutate(Pairwise = hierKernel(z/(w + 0.1)),
         `Higher-order` = kappa*hierKernel(z/(W + 0.1))) |>
  mutate(Total = Pairwise + `Higher-order`) |>
  pivot_longer(!z) |>
  mutate(name = fct_relevel(name, "Total", "Pairwise"),
         type = "With trait variance")

bind_rows(tab1, tab2) |>
  ggplot(aes(x = z, y = value, color = name, alpha = name)) +
  geom_line(linewidth = 0.85) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = c("steelblue", "goldenrod", "forestgreen")) +
  scale_alpha_manual(values = c(0.8, 0.5, 0.5), guide = "none") +
  labs(x = "Trait distance", y = "Competitive effect", color = NULL) +
  facet_wrap(~ type) +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("../../figures/SI/hier_comp_effect.pdf", width = 7, height = 3)
