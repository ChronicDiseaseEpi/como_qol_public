source("Scripts/00_clean_and_transform.R")
rm(cfs)

analysis_main <- "sens" # change between main/sens

if (analysis_main %in% c("main"))  all_df_chsn <- all_df_cmp 
if (analysis_main == "sens")       all_df_chsn <- all_df_cmp_pstv 

## Use SD of EQ5D and SF36 scores to help select consistent priors across outcomes ----
#EQ5D index score for all in UK population data SD=0.23 taken from here https://www.york.ac.uk/che/pdf/DP172.pdf 
#SF36- Based on Norweigan data PCS SD=10.16,  MCS SD=9.08 taken from here https://hqlo.biomedcentral.com/articles/10.1186/s12955-017-0625-9 
# Use same priors as MDM paper for between class, condition and trial variation - https://journals.sagepub.com/doi/10.1177/0272989X211029556
# Use 10 standard deviation for overall effect, this is quite large
my_priors <- c(
  prior(student_t(3, 0, 10), class = Intercept), # population estimate
  prior(student_t(3, 0, 1), class = sd))    # variation of between trial, class and condition prior

## make model formulae 
myforms <- list(
  pool = estimate|se(std.error) ~ 1 + (1|nct_id),
  cond = estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail),
  drug = estimate|se(std.error) ~ 1 + (1|nct_id) + (1|code5),
  dgcn = estimate|se(std.error) ~ 1 + (1|nct_id) + (1|code5) + (1|sail)
)

## EQ5D comorbidity with interaction, original before rescale
if (analysis_main == "main") {
  df_eq5d_nter_orig <- all_df_chsn$df_eq5d_nter
  df_eq5d_nter_orig$estimate <- df_eq5d_nter_orig$estimate*0.23
  df_eq5d_nter_orig$std.error <- df_eq5d_nter_orig$std.error*0.23
  
  nter_brm_drug_cond_orig <- brm(myforms$dgcn, 
                                 data = df_eq5d_nter_orig, 
                                 prior = my_priors, 
                                 save_pars = save_pars(all = TRUE), 
                                 control = list(adapt_delta = 0.99),
                                 cores = 4)
  saveRDS(nter_brm_drug_cond_orig, "nter_explore_orig_scale.Rds")
  
  df_eq5d_como_orig <- all_df_chsn$df_eq5d_como
  df_eq5d_como_orig$estimate <- df_eq5d_como_orig$estimate*0.23
  df_eq5d_como_orig$std.error <- df_eq5d_como_orig$std.error*0.23
  
  brm_drug_cond_orig <- brm(myforms$dgcn, 
                            data = df_eq5d_como_orig, 
                            prior = my_priors, 
                            save_pars = save_pars(all = TRUE), 
                            control = list(adapt_delta = 0.99),
                            cores = 4)
  saveRDS(nter_brm_drug_cond_orig, "explore_orig_scale.Rds")
}

## EQ5D comorbidity with interaction ----
## note that each of these requires recompiling
print(paste0("Scratch_data/", analysis_main, "/", "eq5d_nter_pooled.Rds"))
nter_brm_drug_cond <- brm(myforms$dgcn, data = all_df_chsn$df_eq5d_nter, 
                          prior = my_priors, save_pars = save_pars(all = TRUE), control = list(adapt_delta = 0.99), 
                          cores = 4)
nter_brm_pooled    <- update(nter_brm_drug_cond, formula = myforms$pool)
nter_brm_condition <- update(nter_brm_drug_cond, formula = myforms$cond)
nter_brm_drug      <- update(nter_brm_drug_cond, formula = myforms$drug)

saveRDS(nter_brm_pooled, paste0("Scratch_data/", analysis_main, "/", "eq5d_nter_pooled.Rds"))
saveRDS(nter_brm_condition, paste0("Scratch_data/", analysis_main, "/", "eq5d_nter_cond.Rds"))
saveRDS(nter_brm_drug, paste0("Scratch_data/", analysis_main, "/", "eq5d_nter_drug.Rds"))
saveRDS(nter_brm_drug_cond, paste0("Scratch_data/", analysis_main, "/", "eq5d_nter_drugcond.Rds"))
rm(nter_brm_pooled, nter_brm_condition, nter_brm_drug, nter_brm_drug_cond)

## EQ5D comorbidity without interaction ----
## note that all but "squared" requires recompiling
if (analysis_main == "main") {
  como_brm_drug_cond <- brm(myforms$dgcn, 
                            data = all_df_chsn$df_eq5d_como, 
                            prior = my_priors, 
                            save_pars = save_pars(all = TRUE), 
                            control = list(adapt_delta = 0.99),
                            cores = 4)
  como_sq_eq5d       <- update(como_brm_drug_cond, newdata =  all_df_chsn$df_eq5d_comosq)
  como_brm_pooled    <- update(como_brm_drug_cond, formula = myforms$pool)
  como_brm_condition <- update(como_brm_drug_cond, formula = myforms$cond)
  como_brm_drug      <- update(como_brm_drug_cond, formula = myforms$drug)
  
  saveRDS(como_sq_eq5d, paste0("Scratch_data/", analysis_main, "/", "eq5d_comosq_drugcond.Rds"))
  saveRDS(como_brm_pooled, paste0("Scratch_data/", analysis_main, "/", "eq5d_como_pooled.Rds"))
  saveRDS(como_brm_condition, paste0("Scratch_data/", analysis_main, "/", "eq5d_como_cond.Rds"))
  saveRDS(como_brm_drug, paste0("Scratch_data/", analysis_main, "/", "eq5d_como_drug.Rds"))
  saveRDS(como_brm_drug_cond, paste0("Scratch_data/", analysis_main, "/", "eq5d_como_drugcond.Rds"))
  rm(como_brm_drug_cond, como_sq_eq5d, como_brm_pooled, como_brm_condition, como_brm_drug)
}

## SF36-MCS and PCS all models with drug and condition, note this is the only model which do for como-squared ----  
nter_brm_mcs_drug_cond <- brm(myforms$dgcn, 
                              data = all_df_chsn$df_mcs_nter, 
                              prior = my_priors, 
                              save_pars = save_pars(all = TRUE), 
                              control = list(adapt_delta = 0.999), 
                              cores = 4)
nter_brm_pcs_drug_cond <- update(nter_brm_mcs_drug_cond, newdata = all_df_chsn$df_pcs_nter) # no recompile needed, note, expecting NA statement
saveRDS(nter_brm_mcs_drug_cond, paste0("Scratch_data/", analysis_main, "/", "mcs_nter_drugcond.Rds"))
saveRDS(nter_brm_pcs_drug_cond, paste0("Scratch_data/", analysis_main, "/", "pcs_nter_drugcond.Rds"))

if (analysis_main == "main") {
  como_brm_mcs_drug_cond <- update(nter_brm_mcs_drug_cond, newdata = all_df_chsn$df_mcs_como) # no recompile needed, note, expecting NA statement
  como_brm_pcs_drug_cond <- update(nter_brm_mcs_drug_cond, newdata = all_df_chsn$df_pcs_como) # no recompile needed
  como_sq_mcs            <- update(nter_brm_mcs_drug_cond, newdata = all_df_chsn$df_mcs_comosq) # no recompile needed
  como_sq_pcs            <- update(nter_brm_mcs_drug_cond, newdata = all_df_chsn$df_pcs_comosq, control = list(adapt_delta = 0.999)) # no recompile needed
  saveRDS(como_brm_mcs_drug_cond, paste0("Scratch_data/", analysis_main, "/", "mcs_como_drugcond.Rds"))
  saveRDS(como_brm_pcs_drug_cond, paste0("Scratch_data/", analysis_main, "/", "pcs_como_drugcond.Rds"))
  saveRDS(como_sq_mcs, paste0("Scratch_data/", analysis_main, "/", "mcs_comosq_drugcond.Rds"))
  saveRDS(como_sq_pcs, paste0("Scratch_data/", analysis_main, "/", "pcs_comosq_drugcond.Rds"))
  rm(como_brm_mcs_drug_cond, como_brm_pcs_drug_cond,
     como_sq_mcs, como_sq_pcs)
}
rm(nter_brm_mcs_drug_cond, nter_brm_pcs_drug_cond)

## SF36-MCS and PCS all models with condition  ----
nter_brm_mcs_cond <- brm(myforms$cond, 
                         data = all_df_chsn$df_mcs_nter, 
                         prior = my_priors, 
                         save_pars = save_pars(all = TRUE), 
                         control = list(adapt_delta = 0.999),
                         cores = 4)
nter_brm_pcs_cond <- update(nter_brm_mcs_cond, newdata = all_df_chsn$df_pcs_nter) # no recompile needed, note, expecting NA statement
saveRDS(nter_brm_mcs_cond, paste0("Scratch_data/", analysis_main, "/", "mcs_nter_cond.Rds"))
saveRDS(nter_brm_pcs_cond, paste0("Scratch_data/", analysis_main, "/", "pcs_nter_cond.Rds"))
if (analysis_main == "main") {
  como_brm_mcs_cond <- update(nter_brm_mcs_cond, newdata = all_df_chsn$df_mcs_como) # no recompile needed, note, expecting NA statement
  como_brm_pcs_cond <- update(nter_brm_mcs_cond, newdata = all_df_chsn$df_pcs_como) # no recompile needed
  saveRDS(como_brm_mcs_cond, paste0("Scratch_data/", analysis_main, "/", "mcs_como_cond.Rds"))
  saveRDS(como_brm_pcs_cond, paste0("Scratch_data/", analysis_main, "/", "pcs_como_cond.Rds"))
  rm(como_brm_mcs_cond, como_brm_pcs_cond)
}
rm(nter_brm_mcs_cond, nter_brm_pcs_cond)

## SF36-MCS and PCS all models with drug  ----
nter_brm_mcs_drug <- brm(myforms$drug, 
                         data = all_df_chsn$df_mcs_nter, 
                         prior = my_priors, 
                         save_pars = save_pars(all = TRUE), 
                         control = list(adapt_delta = 0.999),
                         cores = 4)
nter_brm_pcs_drug <- update(nter_brm_mcs_drug, newdata = all_df_chsn$df_pcs_nter) # no recompile needed, note, expecting NA statement
saveRDS(nter_brm_mcs_drug, paste0("Scratch_data/", analysis_main, "/", "mcs_nter_drug.Rds"))
saveRDS(nter_brm_pcs_drug, paste0("Scratch_data/", analysis_main, "/", "pcs_nter_drug.Rds"))
if (analysis_main == "main") {
  como_brm_mcs_drug <- update(nter_brm_mcs_drug, newdata = all_df_chsn$df_mcs_como) # no recompile needed, note, expecting NA statement
  como_brm_pcs_drug <- update(nter_brm_mcs_drug, newdata = all_df_chsn$df_pcs_como) # no recompile needed
  saveRDS(como_brm_mcs_drug, paste0("Scratch_data/", analysis_main, "/", "mcs_como_drug.Rds"))
  saveRDS(como_brm_pcs_drug, paste0("Scratch_data/", analysis_main, "/", "pcs_como_drug.Rds"))
  rm(como_brm_mcs_drug, como_brm_pcs_drug)
}
rm(nter_brm_mcs_drug, nter_brm_pcs_drug)

## SF36-MCS and PCS all models pooled (no drug or condition)  ----
nter_brm_mcs_pool <- brm(myforms$drug, 
                         data = all_df_chsn$df_mcs_nter, 
                         prior = my_priors, 
                         save_pars = save_pars(all = TRUE), 
                         control = list(adapt_delta = 0.999),
                         cores = 4)
nter_brm_pcs_pool <- update(nter_brm_mcs_pool, newdata = all_df_chsn$df_pcs_nter) # no recompile needed, note, expecting NA statement
saveRDS(nter_brm_mcs_pool, paste0("Scratch_data/", analysis_main, "/", "mcs_nter_pooled.Rds"))
saveRDS(nter_brm_pcs_pool, paste0("Scratch_data/", analysis_main, "/", "pcs_nter_pooled.Rds"))
if (analysis_main == "main") {
  como_brm_mcs_pool <- update(nter_brm_mcs_pool, newdata = all_df_chsn$df_mcs_como) # no recompile needed, note, expecting NA statement
  como_brm_pcs_pool <- update(nter_brm_mcs_pool, newdata = all_df_chsn$df_pcs_como) # no recompile needed
  saveRDS(como_brm_mcs_pool, paste0("Scratch_data/", analysis_main, "/", "mcs_como_pooled.Rds"))
  saveRDS(como_brm_pcs_pool, paste0("Scratch_data/", analysis_main, "/", "pcs_como_pooled.Rds"))
  rm(como_brm_mcs_pool, como_brm_pcs_pool)
}
rm(nter_brm_mcs_pool, nter_brm_pcs_pool)
## end----





