# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

library(tidyverse)
library(Hmisc)
library(viridis)
library(patchwork)

library(lme4)
library(MuMIn)

source("00-growth-curve-functions.R")

## Load data ====

eirli <- readRDS("../data/all_data.rds") %>%
  filter(
    dataset == "EIRLI",
    # Those with confirmed no-dx or no follow-up
    (!dx | follow_up == FALSE)
  ) %>%
  group_by(data_id, follow_up, dx) %>%
  select(-dataset, -TOTAL) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow)
  ) %>%
  filter(
    n >= 4
  ) 

eirli_fits <- eirli %>%
  mutate(
    data5 = if_else(n == 5, data, NULL),
    data4 = if_else(n == 4, data, map(data, ~.x[1:4, ])),
    data3 = map(data4, ~.x[1:3, ])
  ) %>%
  select(-data) %>%
  mutate(
    
    max5 = map_if(data5, is.null, ~return(NA), .else = ~max(.x$TOTAL_E)) %>%
      unlist(),
    
    max4 = map_dbl(data4, ~max(.x$TOTAL_E)),
    max3 = map_dbl(data3, ~max(.x$TOTAL_E)),
    
    gc5 = map_if(data5, is.null, ~return(NULL),
                 .else = ~gomp2.fit(.x, response_var = "TOTAL_PLOT", 
                                    t_var = "exact_age", max = 615, 
                                    max.iter = 100)),
    
    gc4 = map(data4, 
              ~gomp2.fit(.x, response_var = "TOTAL_PLOT", t_var = "exact_age", 
                         max = 615, max.iter = 100)),
    gc3 = map(data3, 
              ~gomp2.fit(.x, response_var = "TOTAL_PLOT", t_var = "exact_age", 
                         max = 615, max.iter = 100)),
    
  ) 

# separate this out to avoid running the fits over and over
eirli_fits2 <- eirli_fits %>%
  select(-starts_with("data")) %>%
  pivot_longer(starts_with("gc")) %>%
  filter(
    # NULLs are missing
    !sapply(value, is.null),
    # Length 6 is good data, 1 is a failure to converge
    sapply(value, length) == 6
  ) %>%
  mutate(
    kg = extract.kg(value),
    max = case_when(name == "gc5" ~ max5,
                    name == "gc4" ~ max4,
                    name == "gc3" ~ max3)
  ) %>%
  select(-value, -matches("max[345]")) 

eirli_fits2 %>%
  select(-max) %>%
  pivot_wider(names_from = name, values_from = kg, names_prefix = "kg_") %>%
  ungroup() %>%
  select(starts_with("kg")) %>%
  as.matrix() %>%
  rcorr()

eirli_fits_plot <- eirli_fits2 %>%
  pivot_wider(names_from = name, values_from = kg, names_prefix = "kg_")

# Low group
eirli_fits2 %>%
  select(-max) %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = kg, names_prefix = "kg_") %>%
  filter(
    kg_gc5 < median(kg_gc5, na.rm = TRUE)
  ) %>%
  select(starts_with("kg")) %>%
  as.matrix() %>%
  rcorr()

# High group
eirli_fits2 %>%
  select(-max) %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = kg, names_prefix = "kg_") %>%
  filter(
    kg_gc5 >= median(kg_gc5, na.rm = TRUE)
  ) %>%
  select(starts_with("kg")) %>%
  as.matrix() %>%
  rcorr()

# Difference

eirli_diff <- eirli_fits2 %>%
  mutate(
    name = str_remove(name, "gc")
  ) %>%
  pivot_wider(values_from = c("kg", "max")) %>%
  mutate(
    kg_diff_45 = 100 * (kg_4 - kg_5) / kg_5,
    kg_diff_35 = 100 * (kg_3 - kg_5) / kg_5
  ) %>%
  pivot_longer(starts_with("kg_diff_")) %>%
  na.omit()

png("plots/different_estimates.png", width = 6, height = 4, units = "in", 
       res = 300)

ggplot(eirli_diff, aes(x = max_4, y = value, color = name)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_vline(xintercept = 615 / (exp(1) * 1:2), color = "red") +
  scale_color_manual(values = c("chocolate", "green4"), 
                     labels = c("3-5", "4-5")) +
  scale_x_continuous(breaks = seq(0, 615, by = 100),
    sec.axis = sec_axis(~ 100 * . / 615, name = "% of 615")) +
  theme_bw() +
  labs(x = "Value at 28 months", y = "% error", color = "Contrast") +
  theme(legend.position = "bottom")

dev.off()  

