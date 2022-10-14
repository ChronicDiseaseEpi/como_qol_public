#06_baseline_como_qol
source("Scripts/00_clean_and_transform.R")
cfs_meta <- cfs %>% 
  distinct(nct_id, sail, code5)
rm(all_df, all_df_cmp)

library(tidyverse)

## functions ----
## Can calculate combined standard deviation of groups as per
## https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formula_for_combining_groups.htm
CombMean <- function(n1, m1, n2, m2){
  top <- n1*m1 + n2*m2
  top/(n1 + n2)
}

CombSd <- function(n1, m1, n2, m2, s1, s2){
  top = (n1 - 1)*s1^2 + (n2 - 1)*s2^2 + (n1*n2/(n1 + n2)) * (m1^2 + m2^2 - 2*m1*m2)
  bottom = n1 + n2 - 1
  res = (top/bottom)^0.5
  res
}

CombSdVectorised <- function(n, m, s){
  myrows <- length(n)
  myrows_now <- myrows
  while (myrows_now >= 2) {
    ## select first two values
    n1 <- n[1]
    n2 <- n[2]
    s1 <- s[1]
    s2 <- s[2]
    m1 <- m[1]
    m2 <- m[2]
    ## replace first value with combination of first two values and drop second value
    new_s1 <- CombSd(n1, m1, n2, m2, s1, s2)
    new_m1 <- CombMean(n1, m1, n2, m2)
    new_n1 <- n1 + n2
    s[1] <- new_s1
    m[1] <- new_m1
    n[1] <- new_n1
    s <- s[-2]
    m <- m[-2]
    n <- n[-2]
    # print(s)
    # print(m)
    ## recalculate the length
    myrows_now <- length(s)
  }
  # print(sum(n))
  s
}

## read data ----
base <- read_csv("Data/baseline_summaries_csdr_yoda.csv")

## Note 3 trials had to collapse comorbidity counts as number of participants was less than 5 when exported summary table
collapsed <- c("NCT00231595", "NCT01009086", "NCT01077362")

## Summary statistics ----
## check via simulation that get same mean and SD as via calculations. This is just a check on the coding. The equation comes from the Cochrane handbook ----
recover <- base %>% 
  nest(age = c(age_std_m, age_std_s), outcome = c(base_m, base_s)) 
recover$age     <- map2(recover$participants, recover$age,     ~ rnorm(.x, .y$age_std_m, .y$age_std_s))
recover$outcome <- map2(recover$participants, recover$outcome, ~ rnorm(.x, .y$base_m, .y$base_s))
recover <- recover %>% 
  unnest(cols = c(age, outcome))

recover_smry_age <- recover %>% 
  group_by(como_cnt, parameter) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            age_m = mean(age),
            age_s = sd(age),
            sex = mean(sex),
            outcome_m = mean(outcome),
            outcome_s = sd(outcome)) %>% 
  ungroup() %>% 
  arrange(parameter)
rm(recover)
## calculate using equation from Cochrane handbook ----
base_smry <- base %>% 
  group_by(como_cnt, parameter) %>% 
  ## note can group by any set of parameters
  summarise(age_m = weighted.mean(age_std_m, participants),
            age_s = CombSdVectorised(n = participants, m = age_std_m, s = age_std_s),
            outcome_m = weighted.mean(base_m , participants),
            outcome_s = CombSdVectorised(n = participants, m = base_m, s = base_s),
            male = weighted.mean(sex, participants)) %>%
  ungroup() %>% 
  arrange(parameter) %>% 
  select(como_cnt, parameter, everything())
base_smry_text <- base %>% 
  group_by(nct_id, parameter, arm, sex) %>% 
  summarise(outcome_m = weighted.mean(base_m , participants),
            outcome_s = CombSdVectorised(n = participants, m = base_m, s = base_s),
            n = sum(participants)) %>%
  ungroup()

## By condition
base_smry_class <- base %>% 
  inner_join(cfs_meta) %>% 
  group_by(sail, como_cnt, parameter) %>% 
  ## note can group by any set of parameters
  summarise(trials = sum(!duplicated(nct_id)),
            drug_classes = paste(code5 %>% unique() %>% sort(),collapse = ", "),
            age_m = weighted.mean(age_std_m, participants),
            age_s = CombSdVectorised(n = participants, m = age_std_m, s = age_std_s),
            outcome_m = weighted.mean(base_m , participants),
            outcome_s = CombSdVectorised(n = participants, m = base_m, s = base_s),
            male = paste0(round(100*weighted.mean(sex, participants),1), "%"),
            n = sum(participants)) %>%
  ungroup() %>% 
  arrange(parameter) %>% 
  select(sail, drug_classes, como_cnt, parameter, trials, n, male, everything())

## unstandandardise age 
base_smry_class <- base_smry_class %>% 
  mutate(age_m = round(age_m * 15 + 60),
         age_s = round(age_s*15),
         age = paste0(age_m, " (", age_s, ")"),
         base_m = round(outcome_m, 2),
         base_s = round(outcome_s, 2),
         base = paste0(base_m, " (", base_s, ")")) %>% 
  select(-age_m, -age_s, -outcome_m, -outcome_s,  -base_m, -base_s)
lvls <- c("trials", 
  "n", "age", 
  "male", "base")

base_smry_class_eq5d <- base_smry_class %>% 
  filter(parameter == "eq5d_index") %>% 
  select(-parameter) %>% 
  mutate(como_cnt = paste0("EQ5D, comorbidity count ", como_cnt)) %>% 
  gather("measure", "value", trials:base) %>% 
  pivot_wider(names_from = c(como_cnt), values_from = value, values_fill = "-") %>%
  arrange(sail, drug_classes, factor(measure, lvls))

base_smry_class_sf36 <- base_smry_class %>% 
  filter(parameter == "sf36_pcs") %>% 
  select(-parameter) %>% 
  mutate(como_cnt = paste0("SF36, comorbidity count ", como_cnt)) %>% 
  gather("measure", "value", trials:base) %>% 
  filter(!measure == "base") %>% 
  pivot_wider(names_from = c(como_cnt), values_from = value, values_fill = "-") 
base_smry_class_sf36_bth <- base_smry_class %>% 
  filter(parameter != "eq5d_index") %>% 
  mutate(como_cnt = paste0("SF36, comorbidity count ", como_cnt)) %>% 
  arrange(desc(parameter)) %>% 
  group_by(sail, drug_classes, como_cnt) %>% 
  summarise(value = paste(base, collapse = "; "),
            measure = "base") %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(como_cnt), values_from = value, values_fill = "-") 
base_smry_class_sf36 <- bind_rows(base_smry_class_sf36,
                                  base_smry_class_sf36_bth) %>%
  arrange(sail, drug_classes, factor(measure, lvls))


base_smry_all <- full_join(base_smry_class_eq5d, base_smry_class_sf36) %>% 
  arrange(sail, drug_classes, factor(measure, lvls)) %>% 
  rename(Condition = sail,
         `Drug classes` = drug_classes,
         ` ` = measure) 
write_csv(base_smry_all, "Outputs/baseline_statistics.csv")

RenameFx <- function(x, outcome){
  a <- names(x) %>% 
    str_replace("sail", "Condition") %>% 
    str_replace("drug_classes", "Drug classes")
}

## For paper table 1, need to make more concise ----
## Note 3 trials had to collapse comorbidity counts as number of participants was less than 5 when exported summary table
base_smry_t1 <- base %>% 
  mutate(como_cnt = if_else(nct_id %in% collapsed, NA_integer_, as.integer(como_cnt))) %>% 
  inner_join(cfs_meta) %>% 
  group_by(sail, parameter) %>% 
  ## note can group by any set of parameters
  summarise(trials = sum(!duplicated(nct_id)),
            drug_classes = paste(code5 %>% unique() %>% sort(),collapse = ", "),
            age_m = weighted.mean(age_std_m * 15 + 60, participants),
            age_s = CombSdVectorised(n = participants, m = age_std_m * 15 + 60, s = age_std_s*15),
            outcome_m = weighted.mean(base_m , participants),
            outcome_s = CombSdVectorised(n = participants, m = base_m, s = base_s),
            male = paste0(round(100*weighted.mean(sex, participants),1), "%"),
            n = sum(participants),
            como0 = paste0(round(100*weighted.mean(como_cnt == 0, participants, na.rm = TRUE)), "%"),
            como1 = paste0(round(100*weighted.mean(como_cnt == 1, participants, na.rm = TRUE)), "%"),
            como2 = paste0(round(100*weighted.mean(como_cnt == 2, participants, na.rm = TRUE)), "%")) %>%
  ungroup() %>% 
  arrange(sail, parameter) %>% 
  select(sail, parameter, drug_classes, como0, como1, como2, trials, n, male, everything())

base_smry_t1 <- base_smry_t1 %>% 
  mutate(across(c(age_m:outcome_s), ~ round(.x, 2)),
         outcome = paste0(outcome_m, " (", outcome_s, ")"),
         age = paste0(age_m, " (", age_s, ")")) %>% 
  separate(parameter, into = c("qol", "subqol"), sep = "_") %>% 
  arrange(desc(subqol)) %>% 
  group_by(sail, qol) %>% 
  mutate(
         outcome = paste(outcome, collapse = "; "),
         qol = str_to_upper(qol)) %>% 
  ungroup() %>% 
  filter(!subqol == "mcs") %>% 
  select(sail, qol, drug_classes, trials, n, male, age, como0, como1, como2, `EQ5D or PCS;MCS` = outcome)

write_csv(base_smry_t1 %>% 
            mutate(across(c(como0:como2), ~ if_else(.x == "NaN%", "-", .x))),
          "Outputs/baseline_statistics_paper.csv")

## Summary for figures ----
## By nct_id, sex and arm
base_smry_fig <- base %>% 
  group_by(parameter, nct_id, sex, arm) %>% 
  ## note can group by any set of parameters
  summarise(outcome_m = weighted.mean(base_m , participants),
            outcome_s = CombSdVectorised(n = participants, m = base_m, s = base_s)) %>%
  ungroup()
rm(base, base_smry, base_smry_class, base_smry_class_eq5d, base_smry_class_sf36, base_smry_class_sf36_bth,
   base_smry_all, recover_smry_age)
