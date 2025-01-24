# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:/My Drive/", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

source("00-growth-curve-functions.R")

# Setup ====

library(tidyverse)

wordbank <- readRDS("../../MCDI-analysis/data/Wordbank/WGasWS_WS-scored-230214.rds") %>%
  select(data_id, date_of_test, form, age, comprehension, TOTAL) %>%
  rename(
    production = TOTAL
  ) %>%
  filter(
    !is.na(date_of_test)
  ) %>%
  mutate(
    year = year(date_of_test),
    comprehension = if_else(form == "WS", NA, comprehension),
    year_group = cut(year, 10)
  ) %>%
  pivot_longer(c(comprehension, production)) %>%
  group_by(form, name) %>%
  filter(
    !is.na(value)
  )

wb_nested <- wordbank %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    linear_model = map(data, ~lm(value ~ age + year + age*year, data = .x))
  )

