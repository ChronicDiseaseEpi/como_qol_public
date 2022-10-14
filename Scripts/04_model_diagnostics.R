# Model diagnostics
source("Scripts/02_consolidate_models.R")

## Model diagnostics ----
loo_q <- quietly(loo)
## already have eq5d_inter
loo_res <- map(mdls, ~ loo_q(.x, moment_match = TRUE))
saveRDS(loo_res, "Scratch_data/other_loo.Rds")
loo_res <- readRDS("Scratch_data/other_loo.Rds")

loo_warn <- map(loo_res, ~ tibble(warning = .x$warnings)) %>% bind_rows(.id = "out_term_type") 
loo_rest <- map(loo_res, ~ as_tibble(.x$result$estimates, rownames = "measure")) %>% bind_rows(.id = "out_term_type")
names(loo_rest) <- names(loo_rest ) %>% str_to_lower()

loo_rest <- loo_rest %>% 
  pivot_wider(names_from = measure, values_from = c(estimate, se))
loo_warn <- loo_warn %>% 
  group_by(out_term_type) %>% 
  mutate(sq = seq_along(out_term_type)) %>% 
  ungroup() 
loo_warn <- loo_warn %>% 
  pivot_wider(names_from = sq, values_from = warning, names_prefix = "warn")
loo_rest <- loo_rest %>% 
  left_join(loo_warn)

reloo_slct <- loo_rest %>% 
  filter(!is.na(warn2)) 
reloo_slct$res <- map(reloo_slct$out_term_type, ~ reloo(loo_res[[.x]]$result, mdls[[.x]]))
reloo_slct$rest <- map(reloo_slct$res, ~ as_tibble(.x$estimates, rownames = "measure"))
reloo_slct2 <- reloo_slct %>% 
  select(out_term_type, warn1, warn2, rest) %>% 
  unnest(rest)
names(reloo_slct2) <- names(reloo_slct2) %>% str_to_lower()
reloo_slct2 <- reloo_slct2 %>% 
  pivot_wider(names_from = measure, values_from = c(estimate, se))

loo_final <- bind_rows(loo = loo_rest %>% 
                         filter(is.na(warn2)),
                       reloo = reloo_slct2,
                       .id = "reloo_done")
write_csv(loo_final, "Outputs/loo.csv")
