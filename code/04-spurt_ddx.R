# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

library(tidyverse)

source("00-growth-curve-functions.R")

predict_curve <- function(gc, range = seq(1, 40, 0.001), x) {

  # Create predicted values
  start <- tibble(exact_age = range)
  y_hat <- predict(gc, start)

  end <- start %>%
    mutate(
      # Seems to be some floating point error
      exact_age = round(exact_age, 3),
      y_hat = y_hat
    )

  # Approximate derivative
  steps <- c(1, 0.5, 0.25, 0.1, 0.01, 0.001)

  x_min <- x - steps
  x_max <- x + steps

  y_min <- sapply(steps, function(d) end$y_hat[end$exact_age == round(x - d, 3)],
                  simplify = TRUE)
  y_max <- sapply(steps, function(d) end$y_hat[end$exact_age == round(x + d, 3)],
                  simplify = TRUE)

  rise <- y_max - y_min
  run  <- x_max - x_min

  outcome <- tibble(diff = steps, m = rise / run)

  # Initialize model at starting value (x0) ~ 10, and growth rate (k) = 0
  exp_growth <- nls(m ~ x0 * exp(k * diff), data = outcome,
                    start = list(x0 = 10, k = 0),
                    control = nls.control(minFactor = 10E-6, warnOnly = TRUE))

  # Return limit
  return(unname(coef(exp_growth)["x0"]))

}

good_data <- readRDS("good_data.rds")

good_data2 <- good_data %>%
  unnest(data) %>%
  mutate(
    ddx = map2_dbl(gc, exact_age, ~predict_curve(gc = .x, x = .y))
  )

good_data3 <- good_data2 %>%
  group_by(data_id)  %>%
  mutate(
    max_ddx = max(ddx),
    age_at_max_ddx = exact_age[ddx == max_ddx]
  )

ggplot(good_data3, aes(x = exact_age, y = ddx, color = age_at_max_ddx)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(group = data_id), alpha = 0.1) +
  geom_hline(yintercept = c(10, 21.3), color = "red") +
  geom_vline(xintercept = c(18, 24), color = "blue", linetype = "dashed") +
  scale_color_viridis() +
  facet_wrap(vars(dataset)) +
  theme_bw()

ggplot(good_data3, aes(x = TOTAL_PLOT, y = ddx, color = age_at_max_ddx)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(group = data_id), alpha = 0.1) +
  geom_hline(yintercept = c(10, 21.3), color = "red") +
  scale_color_viridis() +
  facet_wrap(vars(dataset)) +
  theme_bw()

ggplot(filter(good_data3, dataset == "EIRLI"),
       aes(x = TOTAL_PLOT, y = ddx, color = as.factor(age))) +
  geom_point(alpha = 0.2) +
  geom_line(aes(group = data_id), alpha = 0.1) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  geom_hline(yintercept = c(10, 21.3), color = "red") +
  theme_bw()

good_data3 %>%
  group_by(dx_status) %>%
  dplyr::summarize(
    n = n(),
    min_age = min(age_at_max_ddx),
    max_age = max(age_at_max_ddx),
    gt24 = sum(age_at_max_ddx > 24)
  ) %>%
  mutate(
    p_gt24 = round(gt24 / n, 2)
  )

xy <- good_data3 %>%
  mutate(
    inv_at_max_ddx = TOTAL_PLOT[ddx == max_ddx]
  ) %>%
  select(dataset, data_id, dx_status, max_ddx, inv_at_max_ddx, age_at_max_ddx) %>%
  distinct()

library(ggExtra)

(ggplot(xy, aes(x = age_at_max_ddx, y = max_ddx, color = dx_status)) +
  geom_point(alpha = 0.25) +
  geom_smooth(se = FALSE) +
  theme_bw()) %>%
  ggMarginal(type = "density")

(ggplot(xy, aes(x = inv_at_max_ddx, y = max_ddx, color = dx_status)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  theme_bw()) %>%
  ggMarginal(type = "density")
