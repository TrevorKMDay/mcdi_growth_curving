# Setup ====

path <- "/Research/MCDI/growth_curving/code/BEW-aMNLFA"
locs <- c("G:/My Drive", "I:/My Drive", "/Volumes",
          "/home/tkmd/Insync/day00096@umn.edu/Google Drive")

for (i in locs)
  if (dir.exists(paste0(i, path)))
    setwd(paste0(i, path))

# Setup ====

library(tidyverse)
library(aMNLFA)
library(MplusAutomation)

WS_dict <- read_csv("../../data/s_dict.csv")
WS_cat <- unique(WS_dict$category)
WS_dict_count <- WS_dict %>%
  group_by(category) %>%
  summarize(
    n = n()
  )

demo <- readRDS("../../data/all_demographics.rds") %>%
  distinct()


categories <- readRDS("../../data/all_data_by_category.rds") %>%
  select(-starts_with("word_"), -complexity) %>%
  left_join(demo)

length(unique(categories$data_id))


# Analysis ====

eirli_ages <- categories %>%
  filter(
    dataset == "EIRLI"
  ) %>%
  group_by(age) %>%
  dplyr::summarize(
    n = n()
  ) %>%
  filter(
    n > 100,
    age < 30
  ) %>%
  pull(age)

ages_to_compare <- categories %>%
  filter(
    age %in% c(eirli_ages, eirli_ages + 1, eirli_ages - 1)
  ) %>%
  rowwise() %>%
  mutate(
    comparison_age = case_when(
      age %in% eirli_ages ~ age,
      age %in% (eirli_ages - 1) ~ age + 1,
      age %in% (eirli_ages + 1) ~ age - 1
    )
  ) %>%
  pivot_longer(any_of(WS_cat), names_to = "category") %>%
  group_by(dataset, comparison_age, category) %>%
  nest() %>%
  mutate(
    # anova = anova(value ~ dataset),
    n = map_int(data, nrow),
    mean_value = map_dbl(data, ~mean(.x$value)),
    sd_value = map_dbl(data, ~sd(.x$value))
  ) %>% 
  mutate(
    category = str_extract(category, "^...."),
    se_value = sd_value / sqrt(n)
  )

ages_to_compare_test <- categories %>%
  filter(
    age %in% c(eirli_ages, eirli_ages + 1, eirli_ages - 1)
  ) %>%
  rowwise() %>%
  mutate(
    comparison_age = case_when(
      age %in% eirli_ages ~ age,
      age %in% (eirli_ages - 1) ~ age + 1,
      age %in% (eirli_ages + 1) ~ age - 1
    ),
    BCP = dataset == "BCP",
    EIRLI = dataset == "EIRLI"
  ) %>%
  pivot_longer(any_of(WS_cat), names_to = "category") %>%
  group_by(category, comparison_age) %>%
  filter(
    # Missing EIRLI categories
    !(category %in% c("outside", "places", "time_words"))
  ) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(value ~ BCP + EIRLI, data = .x)),
    coef = map(model, ~summary(.x)$coefficients),
    
    bcp_p = map_dbl(coef, ~.x[2, 4]),
    eirli_p = map_dbl(coef, ~.x[3, 4]),
    
    bcp_eff = map_dbl(coef, ~.x[2, 1]),
    eirli_eff = map_dbl(coef, ~.x[3, 1]),
    
    bcp_sig = bcp_p < (.05 / (76 * 2)),
    eirli_sig = eirli_p < (.05 / (76 * 2)),
    
    category = str_extract(category, "^...."),
  )

ages_sig <- ages_to_compare_test %>%
  select(comparison_age, category, ends_with("_sig"), ends_with("_eff")) %>%
  pivot_longer(c(ends_with("_sig"), ends_with("eff"))) %>%
  separate(name, into = c("dataset", "name")) %>%
  mutate(
    dataset = toupper(dataset)
  ) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(
    sig = as.logical(sig),
    eff = if_else(sig, eff, NA_real_)
  )
 
ages_to_compare2 <- ages_to_compare %>%
  left_join(ages_sig) %>%
  filter(
    dataset != "WB",
    !is.na(sig)
  ) 

ggplot(ages_to_compare2, aes(x = interaction(comparison_age, dataset),
                             y = category, fill = eff)) +
  geom_tile() +
  geom_text(aes(label = round(eff, 2))) +
  scale_fill_gradient2() +
  theme_minimal()

# Set up for aMNLFA ====

categories_p <- categories %>%
  pivot_longer(all_of(WS_cat), names_to = "category") %>%
  left_join(WS_dict_count) %>%
  mutate(
    p = value / n
  ) %>%
  pivot_wider(id_cols = c(dataset, inst, data_id, age, exact_age, sex, 
                          mom_college), 
              names_from = category, values_from = p)

cat_aMNLFA <- categories_p  %>%
  rename(
    
    DATASET = dataset,
    INST = inst,
    DATA_ID = data_id,
    
    D_MALE   = sex,
    D_MOMCOL = mom_college,
    
    # ADD SEX
    
    AGE = age,
    AGE_E = exact_age,
    
    L_ACTION = action_words,
    L_ANIMAL = animals,
    L_BODYP  = body_parts,
    L_CLOTH  = clothing,
    L_DESCW  = descriptive_words,
    L_FOODD  = food_drink,
    L_FURNIT = furniture_rooms,
    L_GAMESR = games_routines,
    L_HOUSEH = household,
    L_LOCATE = locations,
    L_OUTSID = outside,
    L_PEOPLE = people,
    L_PLACES = places,
    L_SOUNDS = sounds,
    L_TOYS   = toys,
    L_VEHICL = vehicles,
    
    # Syntax
    S_HELPV  = helping_verbs,
    S_PRON   = pronouns,
    S_QUANT  = quantifiers,
    S_QWORDS = question_words,
    S_TIMEW  = time_words,
    S_CONCTW = connecting_words,
  ) %>%
  mutate(
    # Contrast-code factors
    D_MALE   = if_else(D_MALE == "Male", 1, -1),
    D_MOMCOL = if_else(D_MOMCOL, 1, -1) 
  ) %>%
  select(-L_OUTSID, -L_PLACES, -S_TIMEW, -INST) %>%
  na.omit() %>%
  group_by(DATA_ID) %>%
  mutate(
    # THE ID variable has to be numeric
    ID_NUM = cur_group_id(),
    # Contrast-code factors
    DATA_BCP = if_else(DATASET == "BCP",   1, -1),
    DATA_EIR = if_else(DATASET == "EIRLI", 1, -1)
  ) %>%
  ungroup() %>%
  select(-DATA_ID, -DATASET) 
 
cat_aMNLFA %>%
  filter(
    apply(., 1, function(x) sum(is.na(x))) > 0
  )

one_factor_object <- aMNLFA.object(
  
  # Working directory
  dir = "I:/My Drive/Research/MCDI/growth_curving/code/BEW-aMNLFA/one-factor-dataset/",
  
  # Dataframe
  mrdata = cat_aMNLFA,
  
  # Indicators
  indicators = str_subset(colnames(cat_aMNLFA), "^[LS]_"),
  
  # mean and var are for things you are substantively interested in
  # mean: what your moderators of interest are
  #       Contrast coding of nominal variables
  meanimpact = c("AGE_E", "DATA_BCP", "DATA_EIR", "D_MALE", "D_MOMCOL"),
  
  # var: contrast coding of nominal variables; this is computationally
  #      expensive; JUST DO TIME VARIABLE
  varimpact  = "AGE_E",
  
  # this part: specific indicators impacted by mods? should included all
  #      mean/var impact items
  measinvar  = str_subset(colnames(cat_aMNLFA), "^D"),
  
  factors = c("DATA_BCP", "DATA_EIR", "D_MALE", "D_MOMCOL"),
  
  time = NULL,
  
  ID = "ID_NUM",
  
  # Variables present in DF but not in analysis
  auxiliary = c("AGE"),
  
  # indicate whether you would like to test measurement invariance of
  # thresholds for ordinal indicators. SET TO TRUE. seems to require at
  # least one categorical indicator?
  thresholds = FALSE
  
)

# SETUP =====
set.seed(55455)

## Item plots ====
aMNLFA.itemplots(one_factor_object)

## Calibration sample  ====

file.remove(c("one-factor-dataset/calibration.dat",
              "one-factor-dataset/full.dat"))

aMNLFA.sample(one_factor_object)

# RUN MODEL ====

## Initialize ===

list.files(pattern = "measinvarscript_*") %>%
  file.remove()

# WARNING: varimpactscript.inp has too many characters in the USEVARIABLE line,
# have to edit to fix

aMNLFA.initial(one_factor_object)
runModels("one-factor-dataset/")

aMNLFA.simultaneous(one_factor_object)
runModels("one-factor-dataset/")
