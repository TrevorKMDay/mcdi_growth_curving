# Setup ====

path <- "/Research/MCDI/growth_curving/code"
locs <- c("G:/My Drive", "I:", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

library(tidyverse)

# EIRLI ====

demo_eirli <- readRDS("../data/EIRLI_clean.rds") %>%
  select(data_id, gender, dx, follow_up, exact_age, ends_with("ed")) %>%
    distinct() %>%
    mutate(
      father_college = father_ed %in% c("college", "graduate"),
      mother_college = mother_ed %in% c("college", "graduate"),
    ) %>%
    select(-ends_with("ed"))

demo_eirli_RE <- read_csv("../../MCDI-analysis/data/EIRLI/EIRLI_dx.csv",
                          show_col_types = FALSE) %>%
  select(sid, Ethnicity) %>%
  mutate(
    data_id = sid,
    ethnicity = case_match(Ethnicity, 1 ~ "White", 2 ~ "Black", 3 ~ "API",
                           4 ~ "Hispanic", 5 ~ "AIAN", 6 ~ "Other",
                           NA ~ "Missing")
  ) %>%
  select(data_id, ethnicity) %>%
  left_join(demo_eirli, by = "data_id") %>%
  mutate(
    status = case_when(follow_up & dx ~ "DXplus",
                       follow_up & !dx ~ "DXminus",
                       !follow_up ~ "DX0")
  ) %>%
  select(data_id, status, ethnicity) %>%
  distinct()

demo_eirli_RE %>%
  filter(
    data_id %in% c("E0137", "E0184")
  )

demo_eirli_RE  %>%
  group_by(ethnicity) %>%
  summarize(
    n = n()
  ) %>%
  mutate(
    p = 100 * n / sum(n)
  )

demo_eirli_RE  %>%
  group_by(status, ethnicity) %>%
  summarize(
    n = n()
  ) %>%
  mutate(
    p = signif(100 * n / sum(n), 2)
  ) %>%
  pivot_wider(names_from = status, values_from = c(n, p)) %>%
  janitor::adorn_totals()




demo_eirli %>%
  group_by(follow_up, dx) %>%
  summarize(
    n = n(),
    age_mean = mean(exact_age),
    age_sd = sd(exact_age)
  )

demo_eirli %>%
  select(data_id, follow_up, dx, gender, ends_with("college")) %>%
  distinct() %>%
  group_by(follow_up, dx) %>%
  summarize(
    n = n(),
    male = sum(gender == "Male"),
    fcol = sum(father_college),
    mcol = sum(mother_college)
  ) %>%
  mutate(
    across(c(male, fcol, mcol), ~(.x / n), .names = "{.col}_p")
  )

demo_eirli %>%
  summarize(
    n = n(),
    age_mean = mean(exact_age),
    age_sd = sd(exact_age)
  )

# BCP ====

bcp_data <- readRDS("all_data_curves.rds")  %>%
  filter(dataset == "BCP")

# timepoints
bcp_ages <- bcp_data %>%
  unnest(data)

bcp_ages %>%
  dplyr::summarize(
    m_age = mean(exact_age),
    sd = sd(exact_age)
  )

demo_bcp <- read_csv("../data/demographics-20231213.csv",
                     show_col_types = FALSE) %>%
  select_all(~str_replace(., ".*,", "")) %>%
  filter(
    str_detect(Cohort, "BCP")
  ) %>%
  select(CandID, Sex, ends_with("education"), ends_with("relationship"),
         subject_race, subject_ethnicity) %>%
  mutate(
    data_id = paste0("B", str_pad(CandID, 6, "left", "0"))
  ) %>%
  filter(
    data_id %in% bcp_data$data_id,
    parent1_education != "."
  ) %>%
  distinct() %>%
  mutate(
    subject_ethnicity = replace(subject_ethnicity,
                                subject_ethnicity == "not_answered",
                                NA),
    across(ends_with("education"),
           ~replace(.x, .x == "not_answered", NA) %>%
             str_detect("^college$|grad")),
  )

# write_csv(demo_bcp, "bcp_demo.csv")

demo_bcp <- read_csv("bcp_demo.csv", show_col_types = FALSE)

table(demo_bcp$Sex)

table(demo_bcp$parent1_education, useNA = "a")
table(demo_bcp$parent2_education, useNA = "a")

table(demo_bcp$parent1_relationship, useNA = "a")
table(demo_bcp$parent2_relationship, useNA = "a")

# Wordbank ====

# Load in data
# WS <- read_data("Wordbank/WS-230215.rds")
# WG <- read_data("Wordbank/WG-230215.rds")


# We only need ages for the WS so no need to load WG demographics
WS_demo <- readRDS("../../MCDI-analysis/data/Wordbank/WS-demographics-230215.rds")
WG_demo <- readRDS("../../MCDI-analysis/data/Wordbank/WG-demographics-230215.rds")

WB_demo <- bind_rows(WG_demo, WS_demo) %>%
  select(data_id, form, age, sex, birth_order, caregiver_education, ethnicity, race,
         dataset_name, dataset_name) %>%
  filter(
    !is.na(age)
  )

WB_demo %>%
  group_by(form) %>%
  summarize(
    n = n(),
    age_m = mean(age),
    age_sd = sd(age),

    sex_count = sum(!is.na(sex)),
    male_n = sum(sex == "Male", na.rm = TRUE),

    ed_count = sum(!is.na(caregiver_education)),
    college = sum(caregiver_education %in% c("College", "Some Graduate",
                                             "Graduate"))
  ) %>%
  mutate(
    male_p = male_n / sex_count * 100,
    college_p = college / ed_count * 100,

    sex_p_missing = (n - sex_count) / n * 100,
    ed_p_missing = (n - college) / n * 100
  )

WB_demo %>%
  summarize(
    n = n(),
    age_m = mean(age),
    age_sd = sd(age),

    sex_count = sum(!is.na(sex)),
    male_n = sum(sex == "Male", na.rm = TRUE),

    ed_count = sum(!is.na(caregiver_education)),
    college = sum(caregiver_education %in% c("College", "Some Graduate",
                                             "Graduate"))
  ) %>%
  mutate(
    male_p = male_n / sex_count * 100,
    college_p = college / ed_count * 100,

    sex_p_missing = (n - sex_count) / n * 100,
    ed_p_missing = (n - college) / n * 100
  )

WB_demo %>%
  filter(
    xor(is.na(ethnicity), is.na(race))
  ) %>%
  pull(dataset_name)

WB_RE <- WB_demo %>%
  filter(
    !is.na(race)
  )

table(WB_demo$form,is.na(WB_demo$race))

WB_RE %>%
  group_by(form, ethnicity) %>%
  mutate(
    race = replace_na(as.character(race), "Missing"),
    ethnicity = replace_na(as.character(ethnicity), "Missing")
  ) %>%
  summarize(
    n = n()
  ) %>%
  pivot_wider(names_from = form, values_from = n)

WB_RE_table <- WB_RE %>%
  group_by(form, ethnicity, race) %>%
  mutate(
    race = replace_na(as.character(race), "Missing"),
    ethnicity = replace_na(as.character(ethnicity), "Missing")
  ) %>%
  summarize(
    n = n()
  ) %>%
  pivot_wider(names_from = form, values_from = n) %>%
  mutate(
    total = WG + WS
  )

WB_RE_table %>%
  mutate(
    WG = signif(100 * WG / 3054, 2),
    WS = signif(100 * WS / 5321, 2),
    total = signif(100 * total / (3054 + 5321), 2)
  )


# combine info ====

bcp_demo_combine <- demo_bcp %>%
  select(data_id, Sex, parent1_education) %>%
  rename(
    sex = Sex,
    p_ed = parent1_education
  )

eirli_demo_combine <- demo_eirli %>%
  select(data_id, gender, mother_college) %>%
  rename(
    sex = gender,
    p_ed = mother_college
  ) %>%
  distinct()

wb_demo_combine <- WB_demo %>%
  select(data_id, age,caregiver_education) %>%
  mutate(
    p_ed = caregiver_education %in% c("College", "Some Graduate", "Graduate")
  )

demo_combined <- bind_rows(bcp_demo_combine, eirli_demo_combine,
                           wb_demo_combine)

sex_combined <- table(demo_combined$sex) / sum(!is.na(demo_combined$sex)) * 100
ed_combined <- table(demo_combined$p_ed) / sum(!is.na(demo_combined$p_ed)) * 100

ages <- c(bcp_ages$exact_age, demo_eirli$exact_age, WB_demo$age)

mean(ages)
sd(ages)

## race/ethnicity ====

# White values from ST1.2
white <- c(121 + 1, 1302 + 936, 1640 + 2319, 961)
black <- c(2, 61 + 6 + 125, 69 + 8 + 213, 25)
hispanic <- c(10, 125, 386, 17)

sums <- sapply(list(white, black, hispanic), sum)

# EIRLI mismatch ====

re_ids <- unique(demo_eirli_RE$data_id)
data_ids <- unique(demo_eirli$data_id)

setdiff(re_ids, data_ids)
