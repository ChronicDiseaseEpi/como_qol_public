#01c run unadjusted models

source("Scripts/00_clean_and_transform.R")
rm(cfs)

## load in adjusted models
mdls_names <- c("eq5d_como_cond.Rds", "eq5d_como_drug.Rds", "eq5d_como_drugcond.Rds", 
                "eq5d_como_pooled.Rds", "eq5d_comosq_drugcond.Rds", "eq5d_nter_cond.Rds", 
                "eq5d_nter_drugcond.Rds", "eq5d_nter_drug.Rds", "eq5d_nter_pooled.Rds", 
                "mcs_como_cond.Rds", "mcs_como_drug.Rds", "mcs_como_drugcond.Rds", 
                "mcs_como_pooled.Rds", "mcs_comosq_drugcond.Rds", "mcs_nter_cond.Rds", 
                "mcs_nter_drug.Rds", "mcs_nter_drugcond.Rds", "mcs_nter_pooled.Rds", 
                "pcs_como_cond.Rds", "pcs_como_drug.Rds", "pcs_como_drugcond.Rds", 
                "pcs_como_pooled.Rds", "pcs_comosq_drugcond.Rds", "pcs_nter_cond.Rds", 
                "pcs_nter_drugcond.Rds", "pcs_nter_pooled.Rds", "pcs_nter_drug.Rds")
mdls <- map(mdls_names, ~ readRDS(paste0("Scratch_data/main/", .x)))
mdls_names <- str_sub(mdls_names, 1, -5)
names(mdls) <- mdls_names
str_count(mdls_names, "_") %>% unique()

unad_dfs <- all_df_cmp[str_detect(names(all_df_cmp), "unad$")]
mdls_como <- mdls[str_detect(names(mdls), "como_")]
eq5d <- mdls_como[str_detect(names(mdls_como), "eq5d")]
mcs <- mdls_como[str_detect(names(mdls_como), "mcs")]
pcs <- mdls_como[str_detect(names(mdls_como), "pcs")]
rm(mdls)

# Don't need to recompile these models, same structure, so runs quite quickly
eq5d_unad <- map(eq5d, update, newdata = unad_dfs$df_eq5d_unad)
mcs_unad <- map(mcs, update, newdata = unad_dfs$df_mcs_unad)
pcs_unad <- map(pcs, update, newdata = unad_dfs$df_pcs_unad)
saveRDS(list(eq5d = eq5d_unad,
             mcs = mcs_unad,
             pcs = pcs_unad),
        "Scratch_data/unadjusted_models_como.Rds")

## run models with broader priors on the between trial, between treatment
## comparison, between condition components
## takes a bit longer to run as models need to recompile
my_priors <- c(
  prior(student_t(3, 0, 10), class = Intercept), # population estimate
  prior(student_t(3, 0, 10), class = sd))    # variation of between trial, class and condition prior

mdls_wide <- vector(length = length(mdls_como), mode = "list")
names(mdls_wide) <- names(mdls_como)
for (i in names(mdls_como)) {
  print(i)
  mdls_wide[[i]] <- update(mdls_como[[i]], prior = my_priors)
}
saveRDS(mdls_wide, "Scratch_data/wide_prior.Rds")



effect_smry <- map(mdls_como, function(mdl) posterior_summary(mdl) %>% 
                     as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "out_term_type")
effect_smry_wide <- map(mdls_wide, function(mdl) posterior_summary(mdl) %>% 
                          as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "out_term_type")
effect_smry <- bind_rows(main = effect_smry,
                         wide = effect_smry_wide,
                         .id = "prior_choice")
effect_smry <- effect_smry %>% 
  arrange(out_term_type, params, prior_choice) 
effect_smry <- effect_smry %>% 
  filter(str_detect(params, "code|nct|sail|Intercept"))
write_csv(effect_smry, "Outputs/compare_weak_generic_priors.csv")

cmpr_plot <- ggplot(effect_smry, aes(x = params, 
                                     y = Estimate, ymin = Q2.5, ymax = Q97.5,
                                     colour = prior_choice)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  facet_wrap(~out_term_type, scales = "free_y") +
  theme_clean()
cmpr_plot

effect_smry_diff <- effect_smry %>% 
  group_by(out_term_type, params) %>% 
  summarise(Estimate = (abs(Estimate[1] - Estimate[2]))/abs(Estimate),
          Est.Error =  (abs(Est.Error[1] - Est.Error[2]))/abs(Est.Error)) %>% 
  ungroup()

