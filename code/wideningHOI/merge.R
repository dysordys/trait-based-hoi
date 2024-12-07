library(tidyverse)

clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0) outfile <- clargs[1] else outfile <- "./widen.rds"

tibble(file = Sys.glob("dat_*.rds")) |>
  mutate(dat = map(file, read_rds), .keep = "none") |>
  unnest(dat) |>
  select(!datFile) |>
  write_rds(file = outfile, compress = "xz")
