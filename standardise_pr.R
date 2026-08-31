
library(tidyverse)
library(anatembea)
library(readxl)
library(janitor)
library(epitools)
library(data.table)

## read IPD data

kilifi_trauma_pr <- read_excel("~/Library/CloudStorage/OneDrive-LSTM/KWTRP/Postdoc/KHDSS malaria admission/KCH admissions 1990-2024/data/KCH Trauma PR year_age (160726).xlsx")

### read population data
pop <- read_excel("~/Library/CloudStorage/OneDrive-LSTM/KWTRP/Postdoc/KHDSS malaria admission/KCH admissions 1990-2024/PYO/KHDSS pop 1990-2024_long by age (160726).xlsx")

# ---------------------------------------
# Standard population for the region
# ---------------------------------------

pop <- pop %>% clean_names() %>%
  mutate(yy= year,
         age_yrs = as.numeric(age),
         pop = population) %>%
  select(yy, age_yrs, pop)

# ---------------------------------------
# School-level summaries
# ---------------------------------------

kilifi_trauma_pr2 <- kilifi_trauma_pr %>%
  filter(!is.na(mps_n)) %>%
  group_by(yy,age_yrs) %>%
  summarise(
    sample = n(),
    pos = sum(mps_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rate = pos / sample
  )

# ---------------------------------------
# Merge with standard population
# ---------------------------------------

std_kilifi_trauma_pr <- kilifi_trauma_pr2 %>%
  left_join(pop, by = c("yy", "age_yrs")) %>%
  mutate(
    expected = rate * pop
  )

# ---------------------------------------
# Calculate standardised prevalence
# ---------------------------------------

final_std <- std_kilifi_trauma_pr %>%
  group_by(yy) %>%
  summarise(

    observed_cases = sum(pos, na.rm = TRUE),
    observed_sample = sum(sample, na.rm = TRUE),

    crude_prev = observed_cases / observed_sample * 100,

    expected_cases = sum(expected, na.rm = TRUE),
    total_std_pop = sum(pop, na.rm = TRUE),

    std_prev = expected_cases / total_std_pop * 100,

    .groups = "drop"
  )
