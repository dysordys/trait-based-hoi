library(tidyverse)
library(withr)


randSigmas <- function(S, indVar) {
  if (indVar == "low") {
    runif(S, 0.001, 0.002)
  } else if (indVar == "medium") {
    runif(S, 0.01, 0.03)
  } else {
    runif(S, 0.05, 0.1)
  }
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
crossing(
  replicate = 1:40,
  S = 40,
  w = c(0.1, 0.2, 0.5),
  indVar = c("low", "medium", "high"),
  model = c("evo", "hier", "evoHOI", "hierHOI")
) |>
  rowid_to_column("simID") |>
  mutate(theta = if_else(model %in% c("hier", "hierHOI"), 0, 0.5)) |>
  mutate(sigma = pmap(list(simID, S, indVar),
                      \(id, S, iv) with_seed(9870L*id + 425L, randSigmas(S, iv)))) |>
  mutate(h2 = map2(simID, S,
                   \(id, S) with_seed(12340L*id + 63L, runif(S, 0.1, 0.15)))) |>
  mutate(kappa = case_when(
    model == "evoHOI"  ~ 10,
    model == "hierHOI" ~ 0.2,
    TRUE               ~ NA_real_
  )) |>
  mutate(z0 = if_else(model %in% c("hierHOI"), 0.15, NA_real_)) |>
  mutate(W = if_else(model == "hierHOI", 0.01, NA_real_)) |>
  (\(.) mutate(., pars = pmap(., list)))() |>
  mutate(ic = pmap(list(simID, S, theta, model), \(id, S, theta, model) {
    if (model %in% c("hier", "hierHOI")) {
      c(rep(1, times = S), with_seed(12340L*id + 517L, runif(S, min=0, max=2)))
    } else {
      c(rep(1, times = S), with_seed(12340L*id + 517L, runif(S, min=-theta, max=theta)))
    }
  })) |>
  select(pars, ic, model) |>
  rowid_to_column("datFile") |>
  mutate(datFile = sprintf("dat_%04d.rds", datFile)) |>
  (\(.) split(., .$datFile))() |>
  walk(\(.) write_rds(., .$datFile))

# Submit jobs
tibble(datFile = Sys.glob("dat_*.rds")) |>
  mutate(jobname = str_c("HOI", str_extract(datFile, "\\d+"), ".sh")) |>
  mutate(account = "liu-compute-2023-34", time = "01:30:00", memory = "16000") |>
  mutate(runstr = str_c("simOnCluster.R ", datFile)) |>
  mutate(launch = TRUE) |>
  select(!datFile) |>
  (\(.) pwalk(., bash_script))()
