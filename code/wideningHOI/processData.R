library(tidyverse)


clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 1) {
  inFile <- clargs[1]
  outFile <- clargs[2]
} else {
  inFile <- "../../data/widen.rds"
  outFile <- "../../data/widen-processed.rds"
}


invSimpson <- function(n) {
  p <- n / sum(n)
  1 / sum(p^2)
}


read_rds(inFile) |>
  nest(data = !simID & !model & !W & !replicate) |>
  mutate(diversity = map_dbl(data, \(x) invSimpson(x$n)),
         richness = map_int(data, \(x) length(unique(x$species)))) |>
  select(!data) |>
  write_rds(outFile, compress = "xz")
