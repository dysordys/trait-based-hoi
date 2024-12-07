library(tidyverse)
library(cowplot)



Rcpp::sourceCpp("models-rigged.cpp") # Compile C++ functions implementing the models


hierKernel <- function(x) sapply(x, sigmoid)


gaussKernel <- function(x) exp(-x^2)


# Make a plot of the trait distributions at some given moment
plotSnapshot <- function(sol, bfun = \(z) 1 - exp(-z), sigma = 0.002,
                         limits = c(NA, NA, NA), res = 801) {
  # Normal distribution if the trait z is within 4 sigma of the mean, and NA otherwise:
  truncnorm <- function(z, m, s) if_else(abs(z - m) < 4 * s, dnorm(z, m, s), NA_real_)
  if (is.na(limits[1])) limits[1] <- min(sol$m - 4 * sigma)
  if (is.na(limits[2])) limits[2] <- max(sol$m + 4 * sigma)
  traits <- sol |>
    filter(time == "final state") |>
    crossing(trait = seq(limits[1], limits[2], l = res)) |>
    mutate(density = n * truncnorm(trait, m, sigma)) |>
    mutate(species = as_factor(species))
  ggplot(traits) +
    geom_line(aes(x = trait, y = density, group = species),
              na.rm = TRUE, color = "steelblue") +
    geom_ribbon(aes(x = trait, ymin = 0, ymax = density, group = species),
                fill = "steelblue", alpha = 0.15) +
    geom_function(fun = \(z) max(traits$density, na.rm = TRUE) * bfun(z), n = 501,
                  color = "darkred", alpha = 0.3, na.rm = TRUE, linewidth = 0.75) +
    scale_x_continuous(name = "Trait value", limits = c(limits[1], limits[2])) +
    scale_y_continuous(name = "Density", limits = c(0, limits[3]),
                       expand = expansion(mult = c(0.015, 0.05))) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none",
          plot.margin = unit(c(0.6, 0.3, 0.3, 0.3), "cm"))
}


plotSnaps <- function(data, model) {
  if (str_starts(model, "hier")) {
    plt <- plotSnapshot(data, bfun = \(z) 1 - exp(-z), sigma = 0.004)
  } else {
    plt <- plotSnapshot(data, bfun = \(z) ifelse(abs(z) <= 0.5, 1, 0))
  }
  plt
}


plotKernels <- function(f, w, W, kappa, zmin = -1, zmax = 1, res = 1001) {
  tibble(z = seq(zmin, zmax, l = res)) |>
    mutate(Pairwise = f(z/w), HOI = kappa*f(z/W)) |>
    mutate(Total = Pairwise + HOI) |>
    pivot_longer(!z) |>
    mutate(name = fct_relevel(name, "Total", "Pairwise")) |>
    ggplot(aes(x = z, y = value, color = name, alpha = name)) +
    geom_line(linewidth = 0.85) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = c("steelblue", "goldenrod", "forestgreen")) +
    scale_alpha_manual(values = c(0.8, 0.5, 0.5), guide = "none") +
    labs(x = "Trait distance", y = "Competitive effect", color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}


# Integrate ODE system for a given time sequence
integrateThroughTime <- function(pars, ic, tseq = seq(0, 1e14, l = 11), ...) {
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
  icNew <- tmp$n
  parsNew <- pars
  parsNew$S <- length(survivingSpecies)
  parsNew$m <- pars$m[survivingSpecies]
  parsNew$surv <- survivingSpecies
  final <- integrateToSteady(parsNew, icNew, stol = stol)
  init <- tibble(time = 0, species = 1:pars$S, n = ic, m = pars$m,
                 Sinit = pars$S,w = pars$w, theta = pars$theta,
                 W = pars$W, z0 = pars$z0, kappa = pars$kappa)
  bind_rows(init, final) |>
    filter(n > 1e-5)
}


# Put ODE solution results in tidy tibble
organizeResults <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    rename_with(~str_c(1:pars$S), 1 + 1:pars$S) |>
    pivot_longer(cols = !time, names_to = "species", values_to = "n") |>
    mutate(species = as.integer(species), m = pars$m[species], Sinit = pars$S,
           w = pars$w, theta = pars$theta, W = pars$W, z0 = pars$z0,
           kappa = pars$kappa, replicate = pars$replicate) |>
    mutate(species = pars$surv[species])
}


runModel <- function(model, S, w, W, kappa) {
  tibble(model = model, S = S, w = w, W = W,
         kappa = kappa, theta = 0.5, z0 = 0.15) |>
    mutate(m = map2(S, model, \(S, model) {
      if (str_starts(model, "hier")) seq(0, 2, l = S) else seq(-0.5, 0.5, l = S)
    } )) |>
    mutate(func = case_when(
      model == "evo" ~ list(evo),
      model == "evoHOI" ~ list(evoHOI),
      model == "evoHOI_rigged" ~ list(evoHOI_rigged),
      model == "hier" ~ list(hier),
      model == "hierHOI" ~ list(hierHOI),
      model == "hierHOI_rigged" ~ list(hierHOI_rigged)
    )) |>
    mutate(surv = map(S, \(S) 1:S)) |>
    (\(.) mutate(., pars = pmap(., list)))() |>
    mutate(ic = map(S, \(S) rep(1, times = S))) |>
    select(model, pars, ic) |>
    mutate(sol = map2(pars, ic, solveSteady)) |>
    unnest(sol)
}


runShowModel <- function(model, S, w, W, kappa) {
  runModel(model, S, w, W, kappa) |>
    select(model, time, species, n, m, w, W, z0, kappa) |>
    mutate(time = as_factor(ifelse(time == 0, "initial state", "final state"))) |>
    (\(data) plotSnaps(data, data$model[1]))()
}


runSaveExample <- function(model, w, W, kappa, scen, S = 41) {
  kern <- if (model == "evo") gaussKernel else if (model == "hier") hierKernel else NA
  hoiModel <- str_c(model, "HOI_rigged")
  kernelPlot <- plotKernels(kern, w = w, W = W, kappa = kappa) +
    ggtitle("kernel") + theme(plot.title = element_text(hjust = 0.5))
  pwPlot <- runShowModel(model, S = S, w = w, W = W, kappa = kappa) +
    ggtitle(model) + theme(plot.title = element_text(hjust = 0.5))
  hoiPlot <- runShowModel(hoiModel, S = S, w = w, W = W, kappa = kappa) +
    ggtitle(str_c(model, "HOI")) + theme(plot.title = element_text(hjust = 0.5))
  jointPlot <- plot_grid(pwPlot, kernelPlot, hoiPlot, nrow = 1, align = "hv")
  ggsave(str_c("../../figures/SI/scen-", scen, "-", model, ".pdf"),
         jointPlot, width = 8, height = 3)
}



# Scenario 1: W < w => HOIs do not contribute to diversity:
runSaveExample("evo",  w = 0.3, W = 0.15, kappa = 0.2, scen = 1, S = 41)
runSaveExample("hier", w = 0.3, W = 0.15, kappa = 0.2, scen = 1, S = 41)

# Scenario 2: W > w => HOIs reduce diversity:
runSaveExample("evo",  w = 0.20, W = 0.4, kappa = 0.5, scen = 2, S = 41)
runSaveExample("hier", w = 0.25, W = 0.5, kappa = 5.0, scen = 2, S = 41)

# Scenario 3: W >> 1 => W no longer matters:
runSaveExample("evo",  w = 0.40, W = 1.5, kappa = 0.2, scen = 3, S = 41)
runSaveExample("hier", w = 0.40, W = 3.0, kappa = 0.2, scen = 3, S = 41)

# Scenario 4: W << w but w not too big => extra coexistence:
runSaveExample("evo",  w = 0.40, W = 0.01, kappa = 0.2, scen = 4, S = 41)
runSaveExample("hier", w = 0.80, W = 0.01, kappa = 5.0, scen = 4, S = 41)
