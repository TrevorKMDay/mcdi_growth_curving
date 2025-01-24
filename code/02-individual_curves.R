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
source("00-colors.R")

## Load data ====

site <- read_csv("../data/BCP_Site.csv") %>%
  mutate(
    data_id = paste0("B", str_pad(CandID, 6, "left", "0")),
    site = replace(site, site == "BSB", "UMN")
  ) %>%
  select(-CandID) %>%
  distinct()


all_data <- readRDS("../data/all_data.rds") %>%
  filter(
    dataset != "Wordbank"
  ) %>%
  left_join(site)


all_data_n <- all_data %>%
  group_by(dataset, data_id, follow_up, dx, site) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    dx_status = case_when(
      dataset == "BCP" ~ "BCP",
      !follow_up       ~ "Dx0",
      dx               ~ "Dx+",
      TRUE             ~ "Dx-"
    )
  ) %>%
  filter(
    n >= 3
  ) %>%
  ungroup() %>%
  select(-follow_up, -dx)


length(unique(all_data$data_id[all_data$dataset == "BCP"]))
length(unique(all_data_n$data_id[all_data_n$dataset == "BCP"]))

table(all_data_n$dataset)

table(all_data_n$dx_status, useNA = "a")


# Fit individual curves ====

all_data_curves <- all_data_n %>%
  mutate(
    inv_max = if_else(dataset == "EIRLI", 615, 680),
    data_max = map_dbl(data, ~max(.x$TOTAL_PLOT)),
    max_status = case_when(data_max > inv_max / exp(1) ~ "good",
                           data_max > inv_max / (exp(1) * 2) ~ "backup",
                           TRUE ~ "insufficient"),
    # TOTAL_PLOT is the correct value for the dataset
    gc = map2(data, inv_max,
              ~gomp2.fit(.x, response_var = "TOTAL_PLOT", t_var = "exact_age",
                         max = .y, max.iter = 100)),
    kg = extract.kg(gc)
  )

saveRDS(all_data_curves, "all_data_curves.rds")

all_data_convergence_fail <- all_data_curves %>%
  mutate(
    gc_len = map_int(gc, length)
  ) %>%
  filter(
    gc_len < 6
  ) %>%
  unnest(data)

table(all_data_curves$dataset)

ggplot(all_data_convergence_fail, aes(x = age, y = TOTAL_E)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(vars(data_id)) +
  theme_bw()

## BCP curves ====

bcp_curves <- all_data_curves %>%
  filter(
    dataset == "BCP",
    max_status != "insufficient"
  ) %>%
  select(-data, -dx_status, -inv_max, -gc) %>%
  mutate(
    kU = kg * 680 / exp(1),
    Ti = map_dbl(kg, ~solve.gomp2(y = .x / exp(1), k_g = .x, A = .x + 1)),
    kU_Z = scale(kU)[, 1],
    median_split = if_else(kU >= median(kU), "high", "low")
  )

mean(bcp_curves$Ti)
sd(bcp_curves$Ti)

png("BCP_distribution_site.png", width = 6 , height = 4, units = "in", res = 300)

ggplot(bcp_curves, aes(x = Ti, color = site)) +
  geom_density(linewidth = 1) +
  theme_minimal()

dev.off()

saveRDS(bcp_curves, "bcp_curves.rds")

all_data_curves %>%
  group_by(dataset, n, max_status) %>%
  dplyr::summarize(
    N = n()
  ) %>%
  pivot_wider(names_from = max_status, values_from = N) %>%
  select(dataset, n, insufficient, backup, good) %>%
  write_csv(file = "sample_numbers.csv")

all_data_curves %>%
  filter(
    dataset == "EIRLI"
  ) %>%
  group_by(dx_status, n, max_status) %>%
  dplyr::summarize(
    N = n()
  ) %>%
  pivot_wider(names_from = max_status, values_from = N) %>%
  select(dx_status, n, insufficient, backup, good) %>%
  write_csv(file = "dx_numbers.csv")

# Good data ====

good_data <- all_data_curves %>%
  filter(
    (dataset == "EIRLI" & max_status == "good" & n >= 4) |
      (dataset == "BCP" & max_status != "insufficient")
  ) %>%
  mutate(
    kU = kg * inv_max / exp(1),
    Ti = map2_dbl(inv_max, kg,
                  ~solve.gomp2(y = .x / exp(1), k_g = .y, A = .x + 1)),
    sec_group = if_else(dataset == "BCP", max_status, dx_status),
    site =replace_na(site, "EIRLI"),
    lm = map(data, ~lm(TOTAL_PLOT ~ exact_age, data = .x)),

    gc_AIC = map_dbl(gc, AICc),
    lm_AIC = map_dbl(lm, AICc)

  ) %>%
  na.omit()

ggplot(good_data, aes(x = lm_AIC, y = gc_AIC, color = Ti)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_color_viridis(option = "magma") +
  theme_bw()

 saveRDS(good_data, "good_data.rds")

good_data_n <- good_data  %>%
  group_by(dataset, dx_status, max_status) %>%
  nest()

good_data %>%
  group_by(dataset, dx_status, max_status) %>%
  dplyr::summarize(
    mean_kg = mean(kg),
    mean_kU = mean(kU),
    mean_Ti = mean(Ti)
  )

BCP_ttest <- t.test(good_data_n$data[[1]]$Ti, good_data_n$data[[2]]$Ti)

eirli <- good_data %>%
  filter(
    dataset == "EIRLI"
  ) %>%
  mutate(
    dx0 = dx_status == "Dx0",
    dx_pos = dx_status == "Dx+"
  )

lm(Ti ~ dx0 + dx_pos, data = eirli) %>%
  summary()

png("plots/Ti_values.png", width = 6, height = 3, units = "in", res = 300)

ggplot(good_data, aes(x = dataset, y = Ti)) +
  geom_boxplot(aes(fill = sec_group), outlier.shape = 21,
               outlier.size = 2, outlier.alpha = 0.75) +
  scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, 10),
                     labels = c("Precocious:\n10", "20", "30",
                                "Delayed:\n40")) +
  scale_fill_manual(values = c(colors$BCP$BCP2, colors$EIRLI$DXno,
                               colors$EIRLI$DXyes, colors$EIRLI$DX0,
                               colors$BCP$BCP3),
                    labels = c("Insufficient", "DX-", "DX+", "DX0",
                               "Sufficient"))+
  theme_bw() +
  labs(x = "Dataset", y = bquote(T[i]), fill = "Group", pattern = "Site")

dev.off()

ggplot(good_data, aes(x = dataset, y = kU)) +
  geom_violin(aes(fill = sec_group), draw_quantiles = c(0.25, 0.5, 0.75),
              adjust = 0.75) +
  scale_fill_manual(values = c(colors$BCP$BCP2, colors$EIRLI$DXno,
                               colors$EIRLI$DXyes, colors$EIRLI$DX0,
                               colors$BCP$BCP3),
                    labels = c("Insufficient", "DX-", "DX+", "DX0",
                               "Sufficient"))+
  theme_bw() +
  labs(x = "Dataset", y = bquote(k[U]), fill = "Group", pattern = "Site")

# Individual plots ====

set.seed(55455)

examples <- good_data %>%
  filter(
    sec_group %in% c("good", "Dx-")
  ) %>%
  group_by(dataset) %>%
  slice_sample(n = 2) %>%
  unnest(data)

predict_curve <- function(gc, range) {

  start <- tibble(exact_age = range)
  y_hat <- predict(gc, start)

  end <- start %>%
    mutate(
      y_hat = y_hat
    )

  return(end)

}

example_fits <- good_data %>%
  filter(
    data_id %in% examples$data_id
  ) %>%
  mutate(
    predicted = map(gc, ~predict_curve(.x, 8:40))
  ) %>%
  unnest(predicted)

example_labels <- examples %>%
  select(data_id, kU, Ti) %>%
  distinct() %>%
  mutate(
    across(c(kU, Ti), ~round(.x, 1))
  )

png("./plots/example_participants.png", width = 6, height = 5, units = "in",
    res = 300)

ggplot(examples, aes(x = exact_age, y = TOTAL_PLOT, color = dataset)) +
  geom_point() +
  geom_line(data = example_fits, aes(y = y_hat)) +
  geom_vline(xintercept = 22, color = "blue", size = 1) +
  geom_text(data = example_labels,
            aes(label = paste0("kU=", round(kU, 1))), x = 10, y = 600,
            hjust = 0) +
  geom_text(data = example_labels,
            aes(label = paste0("Ti=", round(Ti, 1))), x = 10, y = 500,
            hjust = 0) +
  facet_wrap(vars(data_id)) +
  scale_x_continuous(limits = c(8, 40)) +
  scale_y_continuous(limits = c(0, 680)) +
  scale_color_manual(values = c(colors$BCP$BCP1, colors$EIRLI$DXno)) +
  theme_bw() +
  labs(x = "Age (mo.)", y = "Inventory total") +
  theme(legend.position = "none")

dev.off()

example_fits %>%
  filter(
    exact_age == 22
  )

# Modeling ====

## EIRLI ====

demo_eirli <- readRDS("../data/EIRLI_clean.rds") %>%
  select(data_id, gender, ends_with("ed")) %>%
  distinct() %>%
  mutate(
    father_college = father_ed %in% c("college", "graduate"),
    mother_college = mother_ed %in% c("college", "graduate"),
  ) %>%
  select(-ends_with("ed"))

good_data_eirli <- good_data %>%
  filter(
    dataset == "EIRLI"
  ) %>%
  left_join(demo_eirli)

eirli_demo_lm <- lm(Ti ~ dx_status + gender + father_college + mother_college,
   data = good_data_eirli)

summary(eirli_demo_lm)$coefficients[, 4] %>%
  p.adjust() %>%
  round(3)

## BCP ====

demo_bcp <- read_csv("../data/BCP-demographics-200609.csv") %>%
  select(CandID, sex, educ_momed, income_inr) %>%
  mutate(
    data_id = paste0("B", str_pad(CandID, 6, "left", "0")),
    mom_college = educ_momed %in% c("college", "grad", "some_grad") %>%
      replace(., . == "not_answered", NA)
  ) %>%
  select(data_id, sex, mom_college, income_inr)

good_data_bcp <- good_data %>%
  filter(
    dataset == "BCP"
  ) %>%
  left_join(demo_bcp)

bcp_demo_lm <- lm(Ti ~ sex + mom_college + income_inr,  data = good_data_bcp)

# BCP site effects ====



good_data_bcp_site <- good_data_bcp %>%
  distinct() %>%
  left_join(site) %>%
  mutate(
    site = replace(site, site == "BSB", "UMN")
  )

ggplot(good_data_bcp_site, aes(x = site, y = kU)) +
  geom_boxplot()

table(good_data_bcp_site$site)

t.test(good_data_bcp_site$kg[good_data_bcp_site$site == "UMN"],
       good_data_bcp_site$kg[good_data_bcp_site$site == "UNC"])

# Plot curves

good_curves <- calculate.curves(good_data$kg, 680, alpha = 0.1)

plot.curves(good_curves, max_val = 680, x_limits = c(0, 40), y_limits = c(0, 680))

# Dissertation figures ====

table(bcp_curves$max_status)

round(mean(bcp_curves$kU), 3)
round(sd(bcp_curves$kU), 3)
round(range(bcp_curves$kU), 3)


round(mean(bcp_curves$Ti), 3)
round(sd(bcp_curves$Ti), 3)
round(range(bcp_curves$Ti), 3)

bcp_curves2 <- left_join(bcp_curves, demo_bcp)

ggplot(bcp_curves2, aes(kU, fill = sex)) +
  geom_boxplot()
