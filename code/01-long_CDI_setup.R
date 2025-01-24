# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:/My Drive/", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

# Setup ====

library(tidyverse)
library(viridis)

source("00-growth-curve-functions.R")

## Load data ====

WS_dict <- read_csv("../data/s_dict.csv")
WS_cat <- unique(WS_dict$category)

eirli_missing_cats <- c("outside", "places", "time_words")

### BCP ====

bcp <- readRDS("../../MCDI-analysis/data/BCP/BCP_all_data_rescored.rds") %>%
  mutate(
    TOTAL_E = select(., any_of(WS_cat), -all_of(eirli_missing_cats)) %>%
      rowSums(),
    dataset = "BCP",
    data_id = paste0("B", data_id)
  ) %>%
  rename(
    exact_age = age,
    age = ideal_age,
  ) %>%
  select(dataset, data_id, inst, age, exact_age, TOTAL, TOTAL_E, everything())


### EIRLI ====

eirli_cats <- readRDS("../data/EIRLI_clean.rds") %>%
  mutate(
    TOTAL = NA,
    TOTAL_E = select(., any_of(WS_cat)) %>%
      rowSums(),
    dataset = "EIRLI",
    inst = "WS"
  ) %>%
  select(dataset, inst, data_id, age, exact_age, follow_up, dx,
         TOTAL, TOTAL_E, everything()) %>%
  rename(
    complexity = COMPLEXITY
  )

eirli <- eirli_cats %>%
  select(dataset, data_id, inst, follow_up, dx, age, exact_age, TOTAL, TOTAL_E)

### Wordbank ====

# Most recent Wordbank
wordbank <- readRDS("../../MCDI-analysis/data/Wordbank/WGasWS_WS-scored-230214.rds") %>%
  filter(
    # Make sure no overlap between WB+EIRLI
    !(dataset_name == "Thal" & form == "WS")
  ) %>%
  mutate(
    dataset = "WB",
    data_id = paste0("WB", data_id),
    # WG as WS in this data
    form = if_else(form == "WG", "WGasWS", "WS"),
    # No more precise data
    exact_age = age,
    TOTAL_E = select(., any_of(WS_cat)) %>%
      rowSums(),
  ) %>%
  rename(
    wb_contributor = dataset_name,
    inst = form
  ) %>%
  select(dataset, data_id, inst, age, exact_age, TOTAL, TOTAL_E, everything())

saveRDS(wordbank, "wordbank_included.rds")

# table(wordbank$dataset_name)

### ALL DATA ====

categories <- bind_rows(bcp, eirli_cats, wordbank) %>%
  select(dataset, data_id, inst, age, exact_age, TOTAL, TOTAL_E,
         all_of(WS_cat))

saveRDS(categories, "../data/all_data_by_category.rds")

all_data <- bind_rows(bcp, wordbank, eirli) %>%
  select(dataset, data_id, inst, age, exact_age, follow_up, dx, TOTAL,
         TOTAL_E) %>%
  mutate(
    TOTAL_PLOT = if_else(dataset == "EIRLI", TOTAL_E, TOTAL)
  )

## Check for drops

all_data2 <- all_data %>%
  arrange(dataset, data_id, exact_age) %>%
  group_by(dataset, data_id, follow_up, dx) %>%
  mutate(
    cum_max = cummax(TOTAL),
    diff_max = cum_max - TOTAL
  )

ids_with_drops <- unique(all_data2$data_id[all_data2$diff_max > 5])

all_data_drops <- all_data2 %>%
  filter(
    data_id %in% ids_with_drops
  ) %>%
  arrange(diff_max)

png("plots/all_data_drops.png", width = 6.5, height = 4,  units = "in", res = 300)

ggplot(all_data_drops, aes(x = exact_age, y = TOTAL, color = inst)) +
  geom_point(aes(shape = exact_age > 18 & inst == "WGasWS"), size = 2.5) +
  geom_line(aes(group = data_id)) +
  facet_wrap(vars(data_id)) +
  theme_bw() +
  labs(x = "Exact age (mo.)", y = "CDI Inventory Total", shape = "Adjusted",
       color = "Instrument")

dev.off()

all_data3 <- all_data2 %>%
  filter(
    # Discard bad individuals
    !(data_id %in% c("B611881", "B999947")),
    # Discard bad timepoint for this person
    !(data_id == "B879509" & age == 30)
  ) %>%
  select(-cum_max, -diff_max) %>%
  ungroup()

rm(all_data, all_data2)
saveRDS(all_data3, "../data/all_data.rds")

all_data3 %>%
  group_by(dataset, inst) %>%
  dplyr::summarize(
    min_exact_age = min(exact_age),
    max_exact_age = max(exact_age)
  )

bcp_no_out_of_range <- bcp %>%
  filter(
    (inst == "WGasWS" & exact_age <= 18) |
      (inst == "WS" & exact_age <= 30)
  )

ggplot(filter(all_data3, exact_age <= 30),
       aes(x = exact_age, y = TOTAL_E, color = dataset)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  geom_smooth(data = bcp_no_out_of_range, color = "black") +
  theme_minimal()

## WG vs WS ====

# Select data occuring in overlapping age range
wg_ws_overlap <- all_data3 %>%
  select(dataset, data_id, inst, age, exact_age, TOTAL) %>%
  filter(
    # EIRLI has weird totals, compare this on only BCP/WB
    dataset != "EIRLI",
    (inst == "WGasWS" & exact_age >= 16) |
      (inst == "WS" & exact_age <= 27.5)
  ) %>%
  arrange(data_id) %>%
  mutate(
    age_c = age - 18,
  )

table(wg_ws_overlap$dataset)

png("plots/wg_ws_overlap.png", width = 5, height = 3,  units = "in",
    res = 300)

ggplot(wg_ws_overlap, aes(x = exact_age, y = TOTAL, color = inst)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(linetype = dataset), linewidth = 1,
              method = "lm", fullrange = FALSE) +
  scale_x_continuous(breaks = 1:30) +
  scale_y_continuous(limits = c(NA, 680), breaks = seq(0, 680, 85)) +
  scale_color_discrete(labels = c("WG", "WS")) +
  coord_cartesian(ylim = c(0, 680)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Age (mo.)", y = "Total words", color = "Instrument")

dev.off()

lm_overlap <- lm(TOTAL ~ age_c + age_c*inst + dataset*inst,
                 data = wg_ws_overlap)
summary(lm_overlap)

lm_overlap_sigterms <- lm(TOTAL ~ age_c + age_c*inst, data = wg_ws_overlap)

all_data4 <- all_data3 %>%
  mutate(
    adj = if_else(inst == "WGasWS",
                  if_else(dataset == "BCP",
                          if_else(exact_age > 18,
                                  coef(lm_overlap_sigterms)[4] +
                                    (exact_age - 18) * sum(coef(lm_overlap_sigterms)[c(2, 4)]),
                                  0),
                          0),
                  0) %>%
      round(),
    TOTAL = TOTAL + adj,
    TOTAL_E = TOTAL_E + round((615 / 680) * adj)
  )

ggplot(filter(all_data4, dataset == "BCP"), aes(x = exact_age, y = TOTAL)) +
  geom_line(aes(group = data_id), alpha = 0.1) +
  geom_point(alpha = 0.25, aes(color = interaction(inst, exact_age > 18)))

rm(bcp, bcp_cats, bcp_demo, bcp_no_out_of_range, bcp_wg, bcp_wg_as_ws, bcp_ws)
rm(all_data_drops, all_data_n, all_demo)
rm(eirli, eirli_cats, eirli_demo)

## Cross-sectional modeling ====

predict_curve <- function(gc, range) {

  start <- tibble(exact_age = range)
  y_hat <- predict(gc, start)

  end <- start %>%
    mutate(
      y_hat = y_hat
    )

  return(end)

}

all_data_n <- all_data4 %>%
  select(-TOTAL_PLOT) %>%
  pivot_longer(starts_with("TOTAL"), names_to = "categories",
               values_to = "TOTAL") %>%
  filter(
    !is.na(TOTAL)
  ) %>%
  group_by(dataset, categories) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    max = if_else(categories == "TOTAL_E", 615, 680),

    # Curve parameters
    gc = map2(data, max,
              ~gomp2.fit(.x, response_var = "TOTAL", t_var = "exact_age",
                         max = .y + 1)),
    y_hat = map(gc, ~predict_curve(.x, 8:40)),

    kg = extract.kg(gc),
    kg_SE = map_dbl(gc, ~summary(.x)$coefficients[1, 2]),
    kg_95CI = 1.96 * kg_SE,

    kU = map2_dbl(kg, max, ~kg * max / exp(1)),
    kU_95CI = kg_95CI * max / exp(1),


    # Calculate time at maximum growth
    Ti = map2_dbl(kg, max, ~solve.gomp2(y = .y / exp(1), k_g = .x, A = .y)),
    Ti_lower = map2_dbl(kg + kg_95CI, max,
                        ~solve.gomp2(y = .y / exp(1), k_g = .x, A = .y)),
    Ti_upper = map2_dbl(kg - kg_95CI, max,
                        ~solve.gomp2(y = .y / exp(1), k_g = .x, A = .y))

  )

all_data_n %>%
  select(dataset, categories, n, kU, kU_95CI, starts_with("Ti")) %>%
  mutate(

    kU_lower = kU - kU_95CI,
    kU_upper = kU + kU_95CI,

    kU_display = paste0(round(kU, 1), " [", round(kU_lower, 1), ", ",
                        round(kU_upper, 1), "]"),
    Ti_display = paste0(round(Ti, 1), " [", round(Ti_lower, 1), ", ",
                        round(Ti_upper, 1), "]")

  )

# x <- predict.nls(object = all_data_n$gc[[1]], newdata = tibble(exact_age = 8:30),
#                 se.fit = TRUE, interval = "confidence", level = 0.95)

predicted_curves <- all_data_n %>%
  select(dataset, categories, y_hat) %>%
  unnest(y_hat)

wordbank_curve <- all_data_n$gc[[4]]

wordbank_predicted <- tibble(exact_age = seq(8, 36, by = 0.001)) %>%
  mutate(
    y_hat = predict(wordbank_curve, newdata = .)
  )

which.min(abs(wordbank_predicted$y_hat - 250))
wordbank_predicted[which.min(abs(wordbank_predicted$y_hat - 250)), ]

ggplot(wordbank_predicted, aes(x = exact_age, y = y_hat)) +
  geom_line()

ggplot(all_data4, aes(x = exact_age, y = TOTAL_E, color = dataset)) +
  geom_point(alpha = 0.1) +
  geom_line(data = filter(predicted_curves, categories == "TOTAL_E"), aes(y = y_hat))

png("plots/all_data.png", width = 5.5, height = 3.3, units = "in",
    res = 300)

ggplot(all_data4, aes(x = exact_age, y = TOTAL_PLOT)) +
  geom_line(aes(group = data_id), alpha = 0.5, color = "grey") +
  geom_point(shape = 21, alpha = 0.2, aes(fill = inst)) +
  geom_line(data = predicted_curves, aes(y = y_hat, linetype = categories),
            color = "black", linewidth = 1) +
  geom_hline(yintercept = c(681, 681 / exp(1)), color = "red") +
  geom_hline(yintercept = c(615, 615 / exp(1)), color = "red",
             linetype = "dashed") +
  scale_x_continuous(limits = c(8, 40), breaks = seq(6, 42, 6)) +
  scale_fill_discrete(labels = c("WG", "WS")) +
  scale_linetype_discrete(labels = c("All", "EIRLI set")) +
  facet_grid(cols = vars(dataset)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Exact age (mo.)", y = "Total words", fill = "Instrument",
       linetype = "Categories")

dev.off()

# Demographics ====

## BCP ====

bcp_demo <- read_csv("../data/BCP-demographics-200609.csv") %>%
  select(CandID, sex, educ_momed) %>%
  rename(
    data_id = CandID
  ) %>%
  mutate(
    data_id = paste0("B", str_pad(data_id, 6, "left", "0")),
    mom_college = educ_momed %in% c("college", "grad", "some_grad") %>%
      replace(., . == "not_answered", NA)
  ) %>%
  select(-educ_momed)

## Wordbank ====

wg_demo <- readRDS("../data/WG-demographics.rds")
ws_demo <- readRDS("../data/WS-demographics.rds")

wb_demo <- bind_rows(wg_demo, ws_demo) %>%
  select(data_id, sex, mom_ed) %>%
  mutate(
    data_id = paste0("WB", data_id),
    mom_college = mom_ed %in% c("College", "Some Graduate", "Graduate")
  ) %>%
  select(-mom_ed)

rm(wg_demo, ws_demo)

## EIRLI ====

eirli_demo <- readRDS("../data/EIRLI_clean.rds") %>%
  select(data_id, gender, mother_ed) %>%
  rename(
    sex = gender,
  ) %>%
  mutate(
    mom_college = mother_ed %in% c("college", "graduate")
  ) %>%
  select(-mother_ed)

all_demo <- bind_rows(bcp_demo, eirli_demo, wb_demo) %>%
  filter(
    data_id %in% all_data3$data_id
  )

saveRDS(all_demo, "../data/all_demographics.rds")

# Grouped analyses

all_data4 <- left_join(all_data3, all_demo)

all_data4_sex_e <- left_join(all_data3, all_demo) %>%
  filter(
    !is.na(sex)
  ) %>%
  group_by(dataset, sex) %>%
  nest() %>%
  mutate(
    gc_E = map(data,
                ~gomp2.fit(.x, response_var = "TOTAL_E", t_var = "exact_age",
                           max = 615)),
    y_hat_e = map(gc_E, ~predict_curve(.x, 8:40)),
  ) %>%
  unnest(y_hat_e)

all_data4_sex_max <- left_join(all_data3, all_demo) %>%
  filter(
    !is.na(sex),
    dataset != "EIRLI"
  ) %>%
  group_by(dataset, sex) %>%
  nest() %>%
  mutate(
    gc = map(data,
               ~gomp2.fit(.x, response_var = "TOTAL", t_var = "exact_age",
                          max = 680)),
    y_hat = map(gc, ~predict_curve(.x, 8:40)),
  ) %>%
  unnest(y_hat)

all_data4_momed_e <- left_join(all_data3, all_demo) %>%
  filter(
    !is.na(mom_college)
  ) %>%
  group_by(dataset, mom_college) %>%
  nest() %>%
  mutate(
    gc = map(data,
               ~gomp2.fit(.x, response_var = "TOTAL_E", t_var = "exact_age",
                          max = 615)),
    y_hat = map(gc, ~predict_curve(.x, 8:40)),
  ) %>%
  unnest(y_hat)

all_data4_momed_max <- left_join(all_data3, all_demo) %>%
  filter(
    !is.na(mom_college),
    dataset != "EIRLI"
  ) %>%
  group_by(dataset, mom_college) %>%
  nest() %>%
  mutate(
    gc = map(data,
             ~gomp2.fit(.x, response_var = "TOTAL", t_var = "exact_age",
                        max = 680)),
    y_hat = map(gc, ~predict_curve(.x, 8:40)),
  ) %>%
  unnest(y_hat)

### Plots ====

sex_e <- ggplot(all_data4_sex_e, aes(x = exact_age, y = y_hat)) +
  geom_line(aes(color = sex, linetype = dataset), linewidth = 1) +
  scale_y_continuous(limits = c(0, 680)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme_bw() +
  labs(x = "Age (mo.)", y = "CDI inventory", color = "Sex",
       linetype = "Dataset") +
  theme(legend.position = "none")

sex_max <- ggplot(all_data4_sex_max, aes(x = exact_age, y = y_hat)) +
  geom_line(aes(color = sex, linetype = dataset), linewidth = 1) +
  scale_y_continuous(limits = c(0, 680)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_bw() +
  labs(x = "Age (mo.)", y = "CDI inventory", color = "Sex",
       linetype = "Dataset") +
  theme(legend.position = "none")

momed_e <- ggplot(all_data4_momed_e, aes(x = exact_age, y = y_hat)) +
  geom_line(aes(color = mom_college, linetype = dataset), linewidth = 1) +
  scale_y_continuous(limits = c(0, 680)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_color_manual(values = c("darkgreen", "darkred")) +
  theme_bw() +
  labs(x = "Age (mo.)", y = "CDI inventory", color = "Mom graduated college",
       linetype = "Dataset") +
  theme(legend.position = "none")

momed_max <- ggplot(all_data4_momed_max, aes(x = exact_age, y = y_hat)) +
  geom_line(aes(color = mom_college, linetype = dataset), linewidth = 1) +
  scale_y_continuous(limits = c(0, 680)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_color_manual(values = c("darkgreen", "darkred")) +
  theme_bw() +
  labs(x = "Age (mo.)", y = "CDI inventory", color = "Mom graduated college",
       linetype = "Dataset") +
  theme(legend.position = "none")

library(patchwork)

(sex_e + sex_max)
