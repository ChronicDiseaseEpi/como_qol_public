# 07_est_interaction_age.R

source("scripts/00_functions_and_packages.R")
library(broom)

## read in data ----
rand <- readRDS("Processed_data/randomisation_for_outcomes.Rds")

## rename rand varaible
rand <- rand %>% 
  rename(subjid = id)

## Load age and comorbidity tables ----
demo <- readRDS(file = "../Extract_Data/Processed_data/all_sponsors_not_conmed.Rds")$demo

## Drop very small number of missing age, values, just exclude these from the data
demo <- demo %>% 
  filter(!is.na(age), !age == "", !age == ".") %>% 
  mutate(age_std = (as.integer(age)-60)/15) %>% 
  select(subjid = id, nct_id, age_std, sex) 

## Add in medical history  and conmed data ----
conmed <- readRDS("../MedDRA_medical_history/Processed_data/medhist_total.Rds")
conmed <- conmed %>% 
  rename(subjid = id) %>% 
  semi_join(rand %>% filter(dm_rndrm_grp %in% c("treat", "control")))

## Calculate agreement ----
conmed_agg <- conmed %>% 
  filter(!is.na(present_conmed), !is.na(present_meddra)) %>% 
  group_by(nct_id, present_conmed, present_meddra, comorbidity) %>% 
  summarise(n = sum(!duplicated(subjid))) %>% 
  ungroup()
meta <- read_csv2("../selected_trials_in_respository_sas_with_medicine_condition.csv")
gsk <- meta %>% 
  filter(sponsor == "GSK") %>% 
  distinct(nct_id) %>% 
  pull(nct_id)

conmed_agg2 <- conmed_agg %>% 
  filter(nct_id %in% gsk) %>% 
  mutate(present_conmed = if_else(present_conmed, "CT", "CF"),
         present_meddra = if_else(present_meddra, "MT", "MF"),
         smry = paste(present_conmed, present_meddra, sep = "|")) %>% 
  select(nct_id, smry, comorbidity, n) %>% 
  spread(smry, n, fill = 0L)
write_csv(conmed_agg2, "Outputs/Agreement_meddra_conmed_definitions.csv")
## Simplify for main anlayis (either TRUE = TRUE) and sensitivity anlaysis (both true = true) ----
## This is deliberately verbose, considering all possibilities
conmed <- conmed %>% 
  mutate(
    main_analysis = case_when(
      is.na(present_conmed) ~ present_meddra,
      is.na(present_meddra) ~ present_conmed,
      present_meddra|present_conmed ~ TRUE,
      (!present_meddra) & (!present_conmed) ~ FALSE
    ),
    sens_analysis = case_when(
      is.na(present_conmed) ~ present_meddra,
      is.na(present_meddra) ~ present_conmed,
      present_meddra & present_conmed ~ TRUE,
      (!present_meddra) | (!present_conmed) ~ FALSE
    )
    ) %>% 
  select(-present_meddra, -present_conmed)

test <- conmed %>% 
  group_by(nct_id, comorbidity) %>% 
  summarise(main = sum(main_analysis),
            sens = sum(sens_analysis)) %>% 
  ungroup()
## appears succesful as can recover numbers with each set
test %>% filter(!main == sens)
test %>% 
  mutate(main_v_sens = case_when(
    is.na(main) & is.na(sens) ~ "All missing",
    main == sens ~ "Only got one",
    main != sens ~ "Got both")) %>% 
  distinct(nct_id, .keep_all = TRUE) %>% 
  count(main_v_sens)

## Drop all missing for analysis
## missing for 14 trials, present for 97
conmed %>% group_by(nct_id) %>% summarise(msng = mean(is.na(main_analysis))) %>% ungroup() %>% count(msng)
conmed <- conmed %>% 
  filter(!is.na(main_analysis))

## Next need to select top 6 conditions for each indication ---- 
top6_cond <- read_csv2("../trials_top6_como.csv") %>% 
  select(nct_id, condition_grp)
top6 <- conmed %>% 
  group_by(nct_id, comorbidity, main_analysis) %>% 
  summarise(n = sum(!duplicated(subjid))) %>% 
  ungroup()
top6 <- top6 %>% 
  mutate(main_analysis = if_else(main_analysis, "x", "n_x")) %>% 
  spread(main_analysis, n, 0L) 
top6_2 <- top6 %>% 
  inner_join(top6_cond) %>% 
  select(-nct_id) %>% 
  group_by(comorbidity, condition_grp) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
  mutate(prop = x/(x + n_x)) %>% 
  arrange(desc(prop)) %>% 
  group_by(condition_grp) %>% 
  slice(1:6) %>% 
  mutate(varname = paste0("v", 1:6)) %>% 
  ungroup() %>% 
  select(condition_grp, comorbidity, varname) %>% 
  inner_join(top6_cond %>% distinct(nct_id, condition_grp)) %>% 
  arrange(condition_grp, nct_id)

## Add in lab and BMI based concomitand data as binaries ----
# lab_cat_lst <- readRDS("../Count_comorbidities/Processed_data/model_covs_categorical.Rds")
# lkp <- c("bmi" = "obese", "dbp" = "CV", "EGFR" = "renal", "fib4" = "liver", "HGB" = "anaem", "sbp" = "CV")
# ## Create set of falses for trials where have the lab data variable
# lab_cat_got <- lab_cat_lst$got_vars
# lab_cat_got <- lab_cat_got %>% 
#   gather("conmed_comp", "got", -trial) %>% 
#   mutate(comorbidity = lkp[conmed_comp]) %>% 
#   filter(got == "Got") %>% 
#   select(trial, comorbidity)
# lab_cat_got <- lab_cat_got %>% 
#   inner_join(rand %>% select(nct_id, trial) %>% distinct())
# rand_got <- rand %>% 
#   distinct(trial, subjid) %>% 
#   inner_join(lab_cat_got) %>% 
#   mutate(lab_analysis = FALSE) %>%
#   select(-trial)
# # Set labe comorbid data in same format
# lab_cat <- lab_cat_lst$all_cov_cat %>% 
#   rename(CV = hyprt) 
# lab_cat <- lab_cat %>% 
#   inner_join(rand %>% rename(id = subjid) %>% select(nct_id, trial, id) %>% distinct()) %>% 
#   select(-trial)
# lab_cat <- lab_cat %>% 
#   gather("comorbidity", "lab_analysis", -nct_id, -id, na.rm = TRUE) %>% 
#   rename(subjid = id)
# ## Drops FALSE where there is a TRUE
# lab_cat <- bind_rows(lab_cat, 
#                      rand_got) %>% 
#   distinct(nct_id, subjid, comorbidity, .keep_all = TRUE)
# rm(lab_cat_lst, rand_got)
# 
# ## Adds CV data from labs. Sets FALSE in conmed/medhist to TRUE if has hypertension
# ## If already TRUE it doesnt change it
# ## It also adds CV FALSE's even if it it just confinded to hypertension (rather than any med)
# conmed_cv <- conmed %>% 
#   filter(comorbidity == "CV", main_analysis == FALSE)
# conmed_cv <- conmed_cv %>% 
#   left_join(lab_cat %>% filter(comorbidity == "CV", lab_analysis == TRUE)) %>% 
#   mutate(main_analysis = if_else(!is.na(lab_analysis), lab_analysis, main_analysis)) %>% 
#   select(-lab_analysis)
# conmed <- conmed %>% 
#   anti_join(conmed_cv %>% select(nct_id, subjid, comorbidity)) %>% 
#   bind_rows(conmed_cv)
# 
# ## Add in non conflincting ones
# conmed <- bind_rows(conmed,
#                      lab_cat %>% filter(!comorbidity == "CV") %>% 
#                        mutate(main_analysis = lab_analysis,
#                               sens_analysis = lab_analysis) %>% 
#                        select(-lab_analysis))
# rm(lab_cat, lab_cat_got)
# 
# ## Dump liver if has an inflammatory condition
meta <- read_csv2("../selected_trials_in_respository_sas_with_medicine_condition.csv") %>%
  distinct(nct_id, condition) %>%
  arrange(condition)
fib4_drp <- c("ankylosing spondylitis", "Axial Spondyloarthritis",
               "Chronic Idiopathic Urticaria (CIU)",
               "Crohn's Disease",  "Psoriasis", "Psoriatic Arthritis",
                "rheumatoid arthritis",  "Ulcerative Colitis and Crihn's disease",
               "Ulcerative Colitis and Crohn's disease")
fib4_drp_nct <- meta$nct_id[meta$condition %in% fib4_drp]
# conmed <- conmed %>% 
#   filter(!( (nct_id %in% fib4_drp_nct) & comorbidity == "liver"))
# 
# ## DUmp if has anaemia, ie same as set to FALSE as this is a count, only loses one
# anaemia <- conmed %>% 
#   filter(comorbidity == "anaem", main_analysis == TRUE) %>%
#   distinct(nct_id, subjid) %>% 
#   mutate(comorbidity = "liver")
# conmed <- conmed %>% 
#   anti_join(anaemia)

## Covert to counts ----
conmed_count_main <- conmed %>% 
  select(-sens_analysis) %>% 
  distinct(nct_id, subjid, comorbidity, .keep_all = TRUE) %>% 
  spread(comorbidity, main_analysis) %>% 
  mutate(skin = FALSE,
         pain = if_else(migraine|arthritis, FALSE, pain),
         arthritis = if_else(inflammatory, FALSE, arthritis))
# liver = if_else(inflammatory|anaem, FALSE, liver)

conmed_count_sens <- conmed %>% 
  select(-main_analysis) %>% 
  distinct(nct_id, subjid, comorbidity, .keep_all = TRUE) %>% 
  spread(comorbidity, sens_analysis) %>% 
  mutate(skin = FALSE,
         pain = if_else(migraine|arthritis, FALSE, pain),
         arthritis = if_else(inflammatory, FALSE, arthritis))
# liver = if_else(inflammatory|anaem, FALSE, liver)

ConmedCount <- function(conmed_count){
  # browser()
  conmed_count <- conmed_count %>% 
    mutate(disease_count = NA) %>%
    ungroup()
  conmed_count$disease_count <- rowSums(conmed_count  %>% select(antacids:urological))  
  conmed_count %>% 
    select(nct_id, subjid, disease_count)
}
conmed_count_main <- ConmedCount(conmed_count_main)
conmed_count_sens <- ConmedCount(conmed_count_sens)
table(conmed_count_main$disease_count, conmed_count_sens$disease_count)

## Add in lab and BMI based as continuous ----
lab_cnt <- readRDS("../Count_comorbidities/Processed_data/model_covs_continuous.Rds")
lab_cnt <- lab_cnt$all_cov_mdl %>% 
  mutate(mbp = sbp*0.5 + dbp*0.5) %>% 
  select(-sbp, -dbp)
lab_cnt <- lab_cnt %>% 
  inner_join(rand %>% rename(id = subjid) %>% select(nct_id, trial, id) %>% distinct()) %>% 
  select(-trial)
lab_cnt %>% 
  arrange(nct_id, id)
bmi_mbp <- lab_cnt %>% 
  select(nct_id, id, bmi, mbp) %>% 
  filter(!is.na(bmi) |!is.na(mbp))
others <- lab_cnt %>% 
  select(nct_id, id, egfr, fib4, hgb) %>% 
  filter(!is.na(egfr) |!is.na(fib4) | !is.na(hgb))
lab_cnt2 <- bmi_mbp %>% 
  inner_join(others)
## select and rename 6 commonest comorbidities
conmed_mdl6 <- conmed %>% 
  inner_join(top6_2)

conmed_mdl6_main <- conmed_mdl6 %>% 
  select(nct_id, subjid, varname, main_analysis) %>% 
  spread(varname, main_analysis)

conmed_mdl6_sens <- conmed_mdl6 %>% 
  select(nct_id, subjid, varname, sens_analysis) %>% 
  spread(varname, sens_analysis)


## Save objects ----
rm(anaemia, conmed, conmed_cv, meta, rand, test, conmed_mdl6)
write_csv(top6_2 %>% spread(varname, comorbidity), "Outputs/six_commonest_conditions.csv")
saveRDS(list(conmed_mdl6_main = conmed_mdl6_main,
             conmed_mdl6_sens = conmed_mdl6_sens,
             conmed_count_main = conmed_count_main,
             conmed_count_sens = conmed_count_sens,
             lab_cnt = lab_cnt2,
             demo = demo),
        "Processed_data/conmed_continuous_categorical.Rds")
