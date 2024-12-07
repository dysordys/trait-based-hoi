library(tidyverse)
library(broom)



testWithCheck <- function(data, method, minDataPoints = 1, ...) {
  if (length(data) > minDataPoints) method(data, ...) else NA
}


compareModels <- function(data, modelTag, property) {
  data |>
    select(model, w, indVar, replicate, {{property}}) |>
    filter(str_starts(model, modelTag)) |>
    pivot_wider(names_from = model, values_from = {{property}}) |>
    mutate(ddiff = get(str_c(modelTag, "HOI")) - get(modelTag), .keep = "unused") |>
    drop_na() |>
    select(!replicate) |>
    nest(ddiff = ddiff) |>
    mutate(ddiff = map(ddiff, unlist)) |>
    mutate(test = map(ddiff, \(x) testWithCheck(x, wilcox.test,
                                                conf.int = FALSE, exact = FALSE))) |>
    mutate(test = map(test, tidy)) |>
    unnest(test) |>
    select(w, indVar, p.value, starts_with("estimate"), starts_with("conf.")) |>
    arrange(w, indVar)
}



dat <- read_rds("../data/data-processed.rds") |>
  mutate(indVar = fct_relevel(indVar, "low", "medium", "high"))

statData <- crossing(tag = c("evo", "hier"),
                     prop = c("diversity", "clustering", "robustness")) |>
  mutate(data = list(dat)) |>
  mutate(tests = pmap(list(data, tag, prop), compareModels)) |>
  unnest(tests)

statData |>
  filter(p.adjust(p.value, "holm") > 0.05) |>
  mutate(tag = str_c(tag, "HOI - ", tag)) |>
  select(!data & !p.value) |>
  rename(`model comparison` = tag,
         `property` = prop,
         `competition width` = w,
         `individual variation` = indVar)
