#06_baseline_qol_model
library(tidyverse)

## read in tranformed data for longitudinal to get class and index condition ----
source("Scripts/00_clean_and_transform.R")
rm(all_df_cmp, cfs, lvl_assign, myoutcome, myoutcome2, myoutcome3)
all_df <- map(all_df, ~ .x %>%  select(nct_id, sail, code5) %>% distinct())
all_df <- all_df %>% 
  bind_rows() %>% 
  distinct()

## Read in baseline models ----
qol_base <- read_csv("Data/baseline_qol_csdr_yoda.csv")

## transform so similar scale, as per longitudinal analysis ----
qol_base <- qol_base %>% 
  mutate(rescale = case_when(
    outcome == "eq5d_index" ~ 0.23,
    outcome == "sf36_mcs" ~ 9.08,
    outcome == "sf36_pcs" ~ 10.16)) %>% 
  mutate(across(c(estimate, se), ~ .x/rescale))

## Add in outcome and class data ----
## Note no increase in rows, nil lost
qol_base <- qol_base %>% 
  inner_join(all_df)

## Select model 4 - age, sex and comorbidity count ----
## checked all makes sense, all have 4 terms per model except two trials which are single-sex trials
qol_base <- qol_base %>% 
  filter(model_n == 4, term == "como_cnt")
eq5d <- qol_base %>% 
  filter(outcome == "eq5d_index")
mcs <- qol_base %>% 
  filter(outcome == "sf36_mcs")
pcs <- qol_base %>% 
  filter(outcome == "sf36_pcs")

## fit brm models ----
my_priors <- c(
  prior(student_t(3, 0, 10), class = Intercept), # population estimate
  prior(student_t(3, 0, 1), class = sd))    # variation of between trial, class and condition prior

eq5d_mdl <- brm(estimate|se(se) ~ 1 + (1|nct_id) + (1|code5) + (1|sail), 
                data = eq5d,
                prior = my_priors, 
                save_pars = save_pars(all = TRUE), 
                control = list(adapt_delta = 0.99))
mcs_mdl <- brm(estimate|se(se) ~ 1 + (1|nct_id) + (1|code5) + (1|sail), 
                data = mcs,
                prior = my_priors, 
                save_pars = save_pars(all = TRUE), 
                control = list(adapt_delta = 0.99))
pcs_mdl <- brm(estimate|se(se) ~ 1 + (1|nct_id) + (1|code5) + (1|sail), 
                data = pcs,
                prior = my_priors, 
                save_pars = save_pars(all = TRUE), 
                control = list(adapt_delta = 0.99))
saveRDS(list(eq5d = eq5d_mdl,
             mcs = mcs_mdl,
             pcs = pcs_mdl),
        "Scratch_data/base_mdls.Rds")
rm(eq5d_mdl, mcs_mdl, pcs_mdl)
mdls <- readRDS("Scratch_data/base_mdls.Rds")
mdls_smry <- map(mdls, ~ summary(.x)$fixed %>% as_tibble())
mdls_smry <- bind_rows(mdls_smry, .id = "outcome")

mdls_smry <- mdls_smry %>% 
  mutate(res = paste0(round(Estimate, 2), " (", round(`l-95% CI`, 2), " to ", round(`u-95% CI`, 2), ")"))
write_csv(mdls_smry, "Outputs/summarise_base_mdls.csv")


## Extract P values WAS ----
ps <- map(mdls, function(x) {
  x <- as_draws_df(x, variable = "b_Intercept")$b_Intercept
  # print(any(p>=0))
  tibble(p = mean(x < 0), n = length(x))
}) %>% 
  bind_rows(.id = "outcome") 
ps <- ps %>% 
  mutate(p = if_else(p == 1, ">0.999", as.character(round(p, 3))))
