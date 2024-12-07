library(tidyverse)
library(withr)



clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 1) {
  inFile <- clargs[1]
  outFile <- clargs[2]
} else {
  inFile <- "../data/data.rds"
  outFile <- "../data/data-processed.rds"
}


invSimpson <- function(n) {
  p <- n / sum(n)
  1 / sum(p^2)
}


traitCV <- function(m) {
  traitdiff <- diff(sort(m))
  sd(traitdiff) / mean(traitdiff)
}


nullCVs <- function(S, replicates) {
  runif(S * replicates) |>
    matrix(replicates, S) |>
    apply(1, traitCV)
}


traitPval <- function(m, replicates) {
  if (length(m) <= 2) { # For less than 3 species, clustering doesn't make sense
    NA_real_
  } else { # For 3 or more species, calculate p-value of clustering:
    mtrans <- (m - min(m)) / (max(m) - min(m)) # Squeeze traits into [0, 1] range
    nulls <- nullCVs(length(mtrans), replicates)
    observed <- traitCV(mtrans)
    length(nulls[nulls < observed]) / replicates
  }
}


sigmoid <- function(x) pnorm(x * sqrt(2)) # Hierarchical competition kernel


epsilonEvoHOI <- function(m, v, w, kappa) {
  S <- length(m)
  eps <- array(0, c(S, S, S))
  for (i in 1:S) for (j in 1:S) for (k in 1:S) {
    dm_ij = (m[i] - m[j])^2
    dm_jk = (m[j] - m[k])^2
    dm_ki = (m[k] - m[i])^2
    num = 8*(dm_jk*v[i] + dm_ki*v[j] + dm_ij*v[k]) + 2*w^2*(dm_ij + dm_jk + dm_ki)
    denom = 16*(v[i]*v[j] + v[j]*v[k] + v[k]*v[i]) + 8*w^2*(v[i] + v[j] + v[k]) + 3*w^4
    eps[i,j,k] <- 2*w*exp(-num/denom) / sqrt(sqrt(pi)*w*denom)
  }
  kappa * eps
}


jac <- function(data, model, w, simID, Sinit = 40) {
  spList <- as.numeric(data$species) # Indices of surviving species
  S <- length(spList)
  z0 <- data$z0[1]
  W <- data$W[1]
  kappa <- data$kappa[1]
  m <- data$m
  dm <- outer(m, m, FUN = `-`)
  ddm <- outer(dm, m, FUN = `-`)
  hddm <- outer(m, outer(m, m, FUN = `+`) / 2, FUN = `-`)
  v <- data$sigma^2
  sv <- outer(v, v, FUN = `+`)
  ssv <- outer(sv, v, FUN = `+`)
  hssv <- outer(v, sv / 4, FUN = `+`)
  if (model %in% c("hier", "hierHOI")) {
    alpha <- sigmoid(dm / sqrt(2*sv + w^2))
  } else {
    alpha <- w*exp(-dm^2 / (2*sv + w^2)) / sqrt(2*sv + w^2)
  }
  if (model == "evoHOI") {
    eps <- epsilonEvoHOI(m, v, w, kappa)
  } else if (model == "hierHOI") {
    eps <- kappa * sigmoid((z0 + hddm) * sqrt(W) / sqrt(2*hssv + W^2))
  } else {
    eps <- array(0, c(S, S, S))
  }
  alphaHOI <- matrix(0, S, S)
  for (i in 1:S) alphaHOI[i,] <- as.numeric((eps[i,,] + t(eps[i,,])) %*% data$n)
  J <- (-data$n) * (alpha + alphaHOI)
  eigen(J, only.values = TRUE)$values |> abs() |> log() |> mean() |> exp()
}



read_rds(inFile) |>
  nest(data = !simID & !model & !indVar & !w & !replicate) |>
  mutate(diversity = map_dbl(data, \(x) invSimpson(x$n)),
         richness = map_int(data, \(x) length(unique(x$species))),
         clustering = map2_dbl(simID, data, \(id, dat)
                               with_seed(49000L*id + 887L, traitPval(dat$m, 1000))),
         robustness = pmap_dbl(list(data, model, w, simID), jac)) |>
  select(!data) |>
  write_rds(outFile, compress = "xz")
