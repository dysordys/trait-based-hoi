library(tidyverse)
library(withr)


randSigmas <- function(S) {
  runif(S, 0.003, 0.009)
}

generate_script <- function(jobname, account, time, memory, runstr) {
  str_c(
    "#!/bin/bash\n",
    "#SBATCH -J ", jobname, "\n",
    "#SBATCH -A ", account, "\n",
    "#SBATCH -t ", time, "\n",
    "#SBATCH --mem=", memory, "\n",
    "#SBATCH -n 1\n",
    "# \n",
    "# Run single task in foreground\n",
    "module add R/4.2.2-hpc1-gcc-11.3.0-bare\n",
    "cd ", Sys.getenv("PWD"), "\n",
    "Rscript ", runstr, "\n",
    "# \n",
    "# Script ends here"
  )
}

bash_script <- function(jobname, account, time, memory, runstr, launch) {
  write(generate_script(jobname, account, time, memory, runstr), jobname)
  if (launch) {
    system(paste("sbatch", jobname))
    Sys.sleep(0.1)
  }
}


# Write temporary input data files
tibble(model = "evoHOI", W = seq(0, 0.11, l = 31)) |>
  bind_rows(tibble(model = "hierHOI", W = seq(0, 1, l = 31))) |>
  crossing(replicate = 1:10, S = 40, w = 0.1) |>
  rowid_to_column("simID") |>
  mutate(theta = if_else(model %in% c("hier", "hierHOI"), 0, 0.5)) |>
  mutate(sigma = map2(simID, S,
                      \(id, S) with_seed(9870L*id + 425L, randSigmas(S)))) |>
  mutate(h2 = map2(simID, S,
                   \(id, S) with_seed(12340L*id + 63L, runif(S, 0.1, 0.15)))) |>
  mutate(m = pmap(list(simID, S, model, theta), \(id, S, model, theta) {
    if (model == "hierHOI") {
      with_seed(12340L*id + 517L, runif(S, min = 0, max = 3))
    } else if (model == "evoHOI") {
      with_seed(12340L*id + 517L, runif(S, min = -theta, max = theta))
    }
  } )) |>
  mutate(kappa = case_when(
    model == "evoHOI"  ~ 10,
    model == "hierHOI" ~ 5,
    TRUE               ~ NA_real_
  )) |>
  mutate(z0 = if_else(model %in% c("hierHOI"), 0.15, NA_real_)) |>
  (\(.) mutate(., pars = pmap(., list)))() |>
  mutate(ic = map(S, \(S) rep(1, times = S))) |>
  select(pars, ic, model) |>
  rowid_to_column("datFile") |>
  mutate(datFile = sprintf("dat_%04d.rds", datFile)) |>
  (\(.) split(., .$datFile))() |>
  walk(\(.) write_rds(., .$datFile))

# Submit jobs
tibble(datFile = Sys.glob("dat_*.rds")) |>
  mutate(jobname = str_c("HOI", str_extract(datFile, "\\d+"), ".sh")) |>
  mutate(account = "liu-compute-2023-34", time = "00:15:00", memory = "16000") |>
  mutate(runstr = str_c("simOnCluster.R ", datFile)) |>
  mutate(launch = TRUE) |>
  select(!datFile) |>
  (\(.) pwalk(., bash_script))()
