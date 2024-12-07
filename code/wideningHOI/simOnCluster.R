library(tidyverse) # Efficient data manipulation and plotting



clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0) datFile <- clargs[1] else datFile <- "dat_0001.rds"


Rcpp::sourceCpp("models-eco.cpp") # Compile C++ functions implementing the models


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
  icNew <- tmp$n
  parsNew <- pars
  parsNew$S <- length(survivingSpecies)
  parsNew$m <- pars$m[survivingSpecies]
  parsNew$sigma <- pars$sigma[survivingSpecies]
  parsNew$h2 <- pars$h2[survivingSpecies]
  integrateToSteady(parsNew, icNew, stol = stol) |> filter(n > 1e-5)
}


# Put ODE solution results in tidy tibble
organizeResults <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    pivot_longer(cols = !time, names_to = "species", values_to = "n") |>
    mutate(species = as_factor(species), m = pars$m[species], Sinit = pars$S, w = pars$w,
           theta = pars$theta, W = pars$W, z0 = pars$z0, kappa = pars$kappa,
           indVar = pars$indVar, sigma = pars$sigma[species], h2 = pars$h2[species],
           replicate = pars$replicate)
}



read_rds(datFile) |>
  mutate(func = case_when(
    model == "evoHOI" ~ list(evoHOI),
    model == "hierHOI" ~ list(hierHOI)
  )) |>
  mutate(pars = map2(pars, func, \(pars, func) { pars$func = func; pars } )) |>
  mutate(simID = map_int(pars, \(x) x$simID), .before = 1) |>
  mutate(sol = map2(pars, ic, solveSteady)) |>
  unnest(sol) |>
  select(!pars & !time) |>
  write_rds(datFile)

cat(str_c("\n\n", datFile, " written successfully\n"))
