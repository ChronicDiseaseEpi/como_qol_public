#02_consolidate_models
## take various model objects saved in "01/run_models.R" and consolidate into a single list

library(tidyverse)
library(brms)
list.files("Scratch_data/main")

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

names(mdls) <- paste0(names(mdls), "_adj")

unad <- readRDS("Scratch_data/unadjusted_models_como.Rds")
unad <- c(unad$eq5d, unad$mcs, unad$pcs)
names(unad) <- paste0(names(unad), "_unad")

mdls <- c(mdls, unad)
rm(unad)

list.files("Scratch_data/sens/")
mdls_names_sens <- c("eq5d_nter_cond.Rds", "eq5d_nter_drug.Rds", "eq5d_nter_drugcond.Rds", 
                     "eq5d_nter_pooled.Rds", "mcs_nter_cond.Rds", "mcs_nter_drug.Rds", 
                     "mcs_nter_drugcond.Rds", "mcs_nter_pooled.Rds", "pcs_nter_cond.Rds", 
                     "pcs_nter_drug.Rds", "pcs_nter_drugcond.Rds", "pcs_nter_pooled.Rds")

mdls_sens <-  map(mdls_names_sens, ~ readRDS(paste0("Scratch_data/sens/", .x)))
mdls_names_sens <- str_sub(mdls_names_sens, 1, -5)
names(mdls_sens) <- mdls_names_sens
