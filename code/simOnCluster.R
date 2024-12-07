library(tidyverse) # Efficient data manipulation and plotting
library(withr) # Setting random number generator locally, without altering it globally



clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0) datFile <- clargs[1] else datFile <- "dat_0001.rds"


Rcpp::sourceCpp("models.cpp") # Compile C++ functions implementing the models


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
                        method = "bdf", stol = 1e-16) {
  tmp <- integrateThroughTime(pars, ic, tseq = tseq, method = method) |>
    filter(time == max(time) & n > 1e-5)
  survivingSpecies <- as.numeric(tmp$species)
  icNew <- c(tmp$n, tmp$m)
  parsNew <- pars
  parsNew$S <- length(survivingSpecies)
  parsNew$sigma <- pars$sigma[survivingSpecies]
  parsNew$h2 <- pars$h2[survivingSpecies]
  integrateToSteady(parsNew, icNew, stol = stol) |> filter(n > 1e-5)
}


# Put ODE solution results in tidy tibble
organizeResults <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    rename_with(~paste0("n_", 1:pars$S), 1 + 1:pars$S) |>
    rename_with(~paste0("m_", 1:pars$S), 1 + pars$S + 1:pars$S) |>
    pivot_longer(cols = !time, names_to = "variable", values_to = "v") |>
    separate(col = variable, into = c("type", "species"), sep = "_", convert = TRUE) |>
    pivot_wider(names_from = "type", values_from = "v") |>
    mutate(species = as_factor(species), Sinit = pars$S, w = pars$w, theta = pars$theta,
           W = pars$W, z0 = pars$z0, kappa = pars$kappa, indVar = pars$indVar,
           sigma = pars$sigma[species], h2 = pars$h2[species],
           replicate = pars$replicate)
}



read_rds(datFile) |>
  mutate(func = case_when(
    model == "evo" ~ list(evo),
    model == "evoHOI" ~ list(evoHOI),
    model == "hier" ~ list(hier),
    model == "hierHOI" ~ list(hierHOI)
  )) |>
  mutate(pars = map2(pars, func, \(pars, func) { pars$func = func; pars } )) |>
  mutate(simID = map_int(pars, \(x) x$simID), .before = 1) |>
  mutate(sol = map2(pars, ic, solveSteady)) |>
  unnest(sol) |>
  select(!pars & !time) |>
  write_rds(datFile)
