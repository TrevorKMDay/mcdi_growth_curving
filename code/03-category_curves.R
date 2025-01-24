# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive",
          )

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

library(tidyverse)

setwd("~/Projects/Minnesota/mcdi_growth_curving/code")

source("00-growth-curve-functions.R")
source("00-colors.R")

## Load data ====

WS_dict <- read_csv("../data/s_dict.csv", show_col_types = FALSE) %>%
  group_by(category) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow)
  )

lexical <- WS_dict$category[1:15]
syntax <- WS_dict$category[16:22]

# Total in each group. For some silly reason sum() and `sum(x) - y` produce
# integers and doubles
lexical_n <- as.numeric(sum(WS_dict$n[1:15]))
lexical_n_E <- lexical_n - 31 - 22
syntax_n <- as.numeric(sum(WS_dict$n[16:22]))
syntax_n_E <- sum(WS_dict$n[16:22])- 12

cats <- readRDS("../data/all_data_by_category.rds")

cats_long <- cats %>%
  pivot_longer(all_of(WS_dict$category), names_to = "category") %>%
  left_join(
    select(WS_dict, category, n)
  ) %>%
  mutate(
    p = value / n,
    inst = if_else(inst == "WS", "WS", "WGasWS"),
    dataset = factor(dataset),
    label = str_to_title(str_replace(category, "_", " "))
  )

png("category_scores_by_study.png", width = 6.5, height = 9, units = "in",
    res = 300)

ggplot(cats_long, aes(x = dataset, y = value, fill = dataset)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE,
              adjust = 1.5) +
  scale_x_discrete(labels = c("BCP", "E", "WB")) +
  scale_y_continuous(n.breaks = 3) +
  coord_flip() +
  facet_wrap(vars(label), scales = "free_x") +
  theme_bw() +
  labs(x = "Dataset", y = "Total words in category known", fill = "Dataset") +
  theme(legend.position = "none", strip.text = element_text(size = 7))

dev.off()

# Category by dataset ====

predict_curve <- function(gc, range) {

  start <- tibble(exact_age = range)
  y_hat <- predict(gc, start)

  end <- start %>%
    mutate(
      y_hat = y_hat
    )

  return(end)

}

cats_long2 <- cats_long %>%
  mutate(
    dx_status = case_when(
      dataset != "EIRLI" ~ NA_character_,
      !follow_up ~ "Dx0",
      dx ~ "Dx+",
      TRUE ~ "Dx-"),
    dataset2 = if_else(is.na(dx_status), dataset, paste0("EIRLI_", dx_status)),
    lexical = category %in% lexical
  )



cats_long2_n <- cats_long2 %>%
  group_by(dataset2, category, n) %>%
  nest() %>%
  mutate(
    gc = map(data,
             ~gomp2.fit(.x, response_var = "p", t_var = "exact_age", max = 1,
                        max.iter = 100)),
    gc_len = map_int(gc, length)
  ) %>%
  filter(
    gc_len == 6
  ) %>%
  mutate(
    pc = map(gc, ~predict_curve(.x, 8:40)),
    kg = extract.kg(gc),
    kU = kg / exp(1),
    kU_words = kU * n,
    Ti = solve.gomp2(y = 1 / exp(1), k_g = kg, A = 1)
  )

cats_stats <- cats_long2_n %>%
  ungroup() %>%
  select(category, dataset2, kU_words, Ti) %>%
  mutate(
    lexical = category %in% lexical
  )

## content/function words ====

content_function_Ti <- cats_stats %>%
  group_by(dataset2, lexical) %>%
  summarize(
    mean_Ti = mean(Ti)
  )

content_function_models <- cats_long2 %>%
  group_by(dataset, dataset2, lexical, data_id, age) %>%
  summarize(
    exact_age = unique(exact_age),
    category_total = sum(value, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(dataset, dataset2, lexical) %>%
  nest() %>%
  mutate(

    # Correct maximum values by group and C/F distinction
    n = if_else(lexical,
                if_else(dataset == "EIRLI", lexical_n_E, lexical_n),
                if_else(dataset == "EIRLI", syntax_n_E, syntax_n)),

    gc = map2(data, n,
             ~gomp2.fit(.x,
                        response_var = "category_total",
                        t_var = "exact_age", max = .y,
                        max.iter = 100)),

    gc_len = map_int(gc, length),

    pc = map(gc, ~predict_curve(.x, 8:40)),

    kg = extract.kg(gc),

    # Calculate kU in units
    kU = kg / exp(1),
    kU_words = kU * n,

    Ti = solve.gomp2(y = 1 / exp(1), k_g = kg, A = 1)
  )

cf_stats <- content_function_models %>%
  ungroup() %>%
  select(dataset2, lexical, kU_words, Ti) %>%
  mutate(
    category = if_else(lexical, "CONTENT", "FUNCTION")
  )

cats_cf_stats <- bind_rows(cats_stats, cf_stats)

# Visualize ====

cats_stats_long <- cats_stats %>%
  pivot_wider(names_from = dataset2, values_from = c(kU_words, Ti)) %>%
  select(category, ends_with("BCP"), ends_with("WB"), ends_with("Dx-"),
         ends_with("Dx0"), ends_with("Dx+"))

write_csv(cats_stats_long, "category_stats_by_dataset.csv")


dir.create("plots/category_growth/", showWarnings = FALSE, recursive = TRUE)

for (i in WS_dict$category) {

  temp <- filter(cats_long, category == i)

  curves <- cats_long2_n %>%
    filter(
      category == i
    )  %>%
    select(-data, -gc) %>%
    unnest(pc)

  p <- ggplot(temp, aes(x = exact_age,  y = p)) +
    geom_point(alpha = 0.25) +
    geom_line(data = curves, aes(y = y_hat, color = dataset2), size = 1.5) +
    geom_hline(yintercept = 1 / exp(1), color = "red", size = 1) +
    facet_wrap(vars(category)) +
    theme_bw()

  png(paste0("plots/category_growth/", i, ".png"), width = 6, height = 4,
      units = "in", res = 300)

  print(p)

  dev.off()

  message(i)

}

category_labels <- c("Sounds", "Animals", "Vehicles", "Toys", "Food/Drink",
                     "Clothing", "Body Parts", "Household Items",
                     "Furniture/Rooms", "Outside Things", "Places to Go",
                     "People", "Games/Routines", "Action Words",
                     "Descriptive Words", "CONTENT WORDS",
                     "Words about Time", "Pronouns",
                     "Question Words",
                     "Prep./Locations", "Quantifiers/Articles",
                     "Helping Verbs", "Connecting Words", "FUNCTION WORDS")

order_for_x_axis <- c(WS_dict$category[1:15], "CONTENT",
                      WS_dict$category[16:22], "FUNCTION")

png("plots/Ti_by_category.png", width = 6, height = 5, units = "in", res = 300)

ggplot(cats_cf_stats, aes(x = category, y = Ti, fill = dataset2)) +
  geom_vline(xintercept = c("CONTENT", "FUNCTION"), linetype = "solid",
             alpha = 0.5, size = 0.75) +
  geom_point(aes(shape = dataset2), size = 3, alpha = 0.75) +
  scale_x_discrete(limits = order_for_x_axis, labels = category_labels) +
  scale_y_continuous(limits = c(14, 33), breaks = seq(14, 33, 2),) +
  scale_fill_manual(values = c(colors$BCP$BCP1, colors$EIRLI$DXno,
                                colors$EIRLI$DXyes, colors$EIRLI$DX0,
                                colors$WB)) +
  scale_shape_manual(values = c(21, 23, 24, 25, 22)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Category", y  = bquote(T[i]), fill = "Group", shape = "Group")

dev.off()

# Individual values ====

lex_v_syn <- cats_long2 %>%
  filter(
    dataset2 != "WB"
  ) %>%
  group_by(dataset2, data_id, age, lexical) %>%
  summarize(
    exact_age = unique(exact_age),
    category_sum = sum(value, na.rm = TRUE)
  ) %>%
  mutate(
    subscore = if_else(lexical, "CONTENT", "FUNCTION")
  ) %>%
  group_by(dataset2, data_id, lexical, subscore) %>%
  nest() %>%
  ungroup() %>%
  select(-lexical) %>%
  mutate(

    n = map_int(data, nrow),
    max_val = map_dbl(data, ~max(.x$category_sum)),

    total = if_else(subscore == "CONTENT",
                    if_else(str_detect(dataset2, "EIRLI"), lexical_n_E, lexical_n),
                    if_else(str_detect(dataset2, "EIRLI"), syntax_n_E, syntax_n)),

    pct_total = max_val / total

  ) %>%
  filter(
    # Sufficient data
    n >= 3,
    pct_total > 1 / (2 * exp(1))
  ) %>%
  mutate(
    gc = map2(data, total,
             ~gomp2.fit(.x,
                          response_var = "category_sum",
                          t_var = "exact_age", max = .y,
                          max.iter = 100)),

    gc_len = map_int(gc, length),
  ) %>%
  filter(
    gc_len > 1
  ) %>%
  mutate(

    pc = map(gc, ~predict_curve(.x, 8:40)),

    kg = extract.kg(gc),

    # Calculate kU in units
    kU = kg / exp(1),
    kU_words = kU * n,

    Ti = solve.gomp2(y = 1 / exp(1), k_g = kg, A = 1)
  )

lex_v_syn_Ti <- lex_v_syn %>%
  select(dataset2, data_id, subscore, Ti) %>%
  pivot_wider(names_from = subscore, values_from = Ti) %>%
  group_by(dataset2) %>%
  nest() %>%
  mutate(
    rcorr = map(data, ~Hmisc::rcorr(as.matrix(.x[, -1]))),
    r = map_dbl(rcorr, ~.x$r[1, 2]),
    n = map_dbl(rcorr, ~.x$n[1, 2]),
    p = map_dbl(rcorr, ~.x$P[1, 2]),
  )

cor(lex_v_syn_Ti$CONTENT, lex_v_syn_Ti$FUNCTION, use = "c")
