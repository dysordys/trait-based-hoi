library(tidyverse)
library(withr)



Rcpp::sourceCpp("models.cpp") # Compile C++ functions implementing the models


# Make a plot of the trait distributions at some given moment
plotSnapshot <- function(sol, bfun = \(z) 1 - exp(-z), limits = c(NA,NA,NA), res = 501) {
  # Normal distribution if the trait z is within 4 sigma of the mean, and NA otherwise:
  truncnorm <- function(z, m, s) if_else(abs(z - m) < 4 * s, dnorm(z, m, s), NA_real_)
  if (is.na(limits[1])) limits[1] <- min(sol$m - 4 * sol$sigma)
  if (is.na(limits[2])) limits[2] <- max(sol$m + 4 * sol$sigma)
  traitaxis <- seq(limits[1], limits[2], l = res)
  traits <- crossing(sol, trait = traitaxis) |>
    mutate(density = n * truncnorm(trait, m, sigma))
  landscape <- tibble(trait = traitaxis) |> # for plotting intrinsic rates
    mutate(r = bfun(trait) * if_else(
      is.na(limits[3]), max(traits$density, na.rm = TRUE), limits[3] * 0.98)) |>
    mutate(r = if_else(r < 0, NA_real_, r)) # where b(z) < 0: set to NA
  traits |>
    mutate(species = as_factor(species)) |>
    ggplot() +
    geom_line(aes(x = trait, y = density, colour = species), na.rm = TRUE) +
    geom_ribbon(aes(x = trait, ymin = 0, ymax = density, fill = species), alpha = 0.15) +
    geom_line(data = landscape, aes(x = trait, y = r),
              colour = "darkred", alpha = 0.5, na.rm = TRUE, linewidth = 0.75) +
    scale_x_continuous(name = "trait value", limits = c(limits[1], limits[2])) +
    scale_y_continuous(name = "density", limits = c(0, limits[3]),
                       expand = expansion(mult = c(0.015, 0.05))) +
    facet_wrap(~ time, ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(0.6, 0.3, 0.3, 0.3), "cm"))
}


plotSnaps <- function(data, model, limits = c(NA,NA,NA), res = 501) {
  if (model %in% c("hier", "hierHOI")) {
    plt <- plotSnapshot(data, bfun = \(z) 1 - exp(-z), limits, res)
  } else {
    plt <- plotSnapshot(data, bfun = \(z) ifelse(abs(z) <= 0.5, 1, 0), limits, res)
  }
  if (model %in% c("evo", "evoHOI")) plt else plt + theme(axis.title.y = element_blank())
}


# Integrate ODE system for a given time sequence
integrateThroughTime <- function(pars, ic, tseq, ...) {
  deSolve::ode(func = pars$func, y = ic, parms = pars, times = tseq, ...) |>
    organizeResults(pars)
}


# Integrate ODE system until its fixed point is reached
integrateToSteady <- function(pars, ic, ...) {
  rootSolve::runsteady(y = ic, time = c(0, Inf), func = pars$func,
                       parms = pars, ...)$y |>
    enframe() |>
    pivot_wider() |>
    mutate(time = Inf, .before = 1) |>
    organizeResults(pars)
}


solveSteady <- function(pars, ic, tseq = seq(0, 1e14, l = 11),
                        method = "bdf", stol = 1e-15) {
  tmp <- integrateThroughTime(pars, ic, tseq = tseq, method = method) |>
    filter(time == max(time) & n > 1e-5)
  survivingSpecies <- as.numeric(tmp$species)
  icNew <- c(tmp$n, tmp$m)
  parsNew <- pars
  parsNew$S <- length(survivingSpecies)
  parsNew$sigma <- pars$sigma[survivingSpecies]
  parsNew$h2 <- pars$h2[survivingSpecies]
  final <- integrateToSteady(parsNew, icNew, stol = stol)
  init <- tibble(
    time = 0,
    species = 1:pars$S,
    n = ic[1:pars$S],
    m = ic[(pars$S + 1):length(ic)],
    Sinit = pars$S, w = pars$w, theta = pars$theta, W = pars$W, z0 = pars$z0,
    kappa = pars$kappa, sigma = pars$sigma, h2 = pars$h2
  )
  bind_rows(init, final) |>
    filter(n > 1e-5)
}


# Put ODE solution results in tidy tibble
organizeResults <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    rename_with(~paste0("n_", 1:pars$S), 1 + 1:pars$S) |>
    rename_with(~paste0("m_", 1:pars$S), 1 + pars$S + 1:pars$S) |>
    pivot_longer(cols = !time, names_to = "variable", values_to = "v") |>
    separate(col = variable, into = c("type", "species"), sep = "_", convert = TRUE) |>
    pivot_wider(names_from = "type", values_from = "v") |>
    mutate(species = species, Sinit = pars$S, w = pars$w, theta = pars$theta,
           W = pars$W, z0 = pars$z0, kappa = pars$kappa, indVar = pars$indVar,
           sigma = pars$sigma[species], h2 = pars$h2[species],
           replicate = pars$replicate)
}



tibble(
  model = c("evo", "hier", "evoHOI", "hierHOI"),
  S = 40,
  w = 0.1
) |>
  rowid_to_column("simID") |>
  mutate(theta = if_else(model %in% c("hier", "hierHOI"), 0, 0.5)) |>
  mutate(sigma = map2(simID, S,
                      \(id, S) with_seed(6870L*id + 425L, runif(S, 0.003, 0.009)))) |>
  mutate(h2 = map2(simID, S,
                   \(id, S) with_seed(14340L*id + 63L, runif(S, 0.1, 0.15)))) |>
  mutate(kappa = case_when(
    model == "evoHOI"  ~ 10,
    model == "hierHOI" ~ 5,
    TRUE               ~ NA_real_
  )) |>
  mutate(z0 = if_else(model %in% c("hierHOI"), 0.15, NA_real_)) |>
  mutate(W = if_else(model == "hierHOI", 0.01, NA_real_)) |>
  (\(.) mutate(., pars = pmap(., list)))() |>
  mutate(ic = pmap(list(simID, S, theta, model), \(id, S, theta, model) {
    init_n <- case_when(
      model == "evo"     ~ rep(1, times = S),
      model == "hier"    ~ rep(0.1, times = S),
      model == "evoHOI"  ~ rep(0.1, times = S),
      model == "hierHOI" ~ rep(0.1, times = S),
    )
    if (model %in% c("hier", "hierHOI")) {
      c(init_n, with_seed(14340L*id + 517L, runif(S, min = 0, max = 2)))
    } else {
      c(init_n, with_seed(14340L*id + 517L, runif(S, min = -theta, max = theta)))
    }
  })) |>
  select(simID, model, pars, ic) |>
  mutate(func = case_when(
    model == "evo" ~ list(evo),
    model == "evoHOI" ~ list(evoHOI),
    model == "hier" ~ list(hier),
    model == "hierHOI" ~ list(hierHOI)
  )) |>
  mutate(pars = map2(pars, func, \(pars, func) { pars$func = func; pars } )) |>
  mutate(sol = map2(pars, ic, solveSteady)) |>
  unnest(sol) |>
  select(model, time, species, n, m, Sinit, w, W, z0, kappa, sigma, h2) |>
  write_rds("../data/examples.rds", compress = "xz")


read_rds("../data/examples.rds") |>
  mutate(time = as_factor(ifelse(time == 0, "initial state", "final state"))) |>
  nest(data = !model) |>
  mutate(model = fct_relevel(model, "evo", "hier", "evoHOI", "hierHOI")) |>
  arrange(model) |>
  mutate(plot = map2(data, model, \(d, m) plotSnaps(d, m, res = 801))) |>
  pull(plot) |>
  (\(.) cowplot::plot_grid(plotlist = ., nrow = 2, align = "hv", vjust = 1.0,
                           labels = c("evo", "hier", "evoHOI", "hierHOI")))()
ggsave("../figures/examples.pdf", width = 6, height = 6)
