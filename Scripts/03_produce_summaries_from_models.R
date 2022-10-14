# 03 produce summaries from model
source("Scripts/00_clean_and_transform.R")
# do not run ("Scripts/01_run_models.R") each time as time consuming
source("Scripts/02_consolidate_models.R")
library(tidybayes)
library(ggthemes)

## Extract estimates for overall effects and standard deviations ----
effect_smry <- map(mdls, function(mdl) posterior_summary(mdl) %>% 
                     as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "out_term_type")
effect_smry <- effect_smry %>% 
  filter(params %in% c("b_Intercept",  "sd_code5__Intercept", "sd_nct_id__Intercept", 
                       "sd_sail__Intercept"))


## Run this function so we can save out warnings into a dataframe
summary_q <- quietly(summary)
warn_smry <- map(mdls, ~ tibble(warning = summary_q(.x)$warnings)) %>% 
  bind_rows(.id = "out_term_type")
effect_smry <- effect_smry %>% 
  left_join(warn_smry) %>% 
  SepModelNames()
## add code to cope with adjusted and unadjusted

write_csv(effect_smry, "Outputs/summarise_long.csv")

t2 <- effect_smry %>% 
  filter(params == "b_Intercept", term == "como") %>% 
  select(outcome, model_type, adj, Estimate,  Q2.5,   Q97.5) %>% 
  mutate(across(c(Estimate,  Q2.5,   Q97.5), ~ round(.x, 2))) %>% 
  mutate(res = paste0(Estimate, " (", Q2.5, " to ", Q97.5, ")")) %>% 
  select(outcome, model_type, adj, res)

t3 <- effect_smry %>% 
  filter(params == "b_Intercept", term == "nter") %>% 
  select(outcome, model_type, adj, Estimate,  Q2.5,   Q97.5) %>% 
  mutate(across(c(Estimate,  Q2.5,   Q97.5), ~ round(.x, 4))) %>% 
  mutate(res = paste0(Estimate, " (", Q2.5, " to ", Q97.5, ")")) %>% 
  select(outcome, model_type, adj, res) %>% 
  spread(outcome, res)

t2_unscale <- effect_smry %>% 
  filter(params == "b_Intercept", term == "como") %>% 
  select(outcome, model_type, adj, Estimate,  Q2.5,   Q97.5) %>% 
  mutate(across(c(Estimate,  Q2.5,   Q97.5), ~ case_when(
    outcome == "eq5d" ~ .x*0.23,
    outcome == "mcs" ~ .x*9.08,
    outcome == "pcs" ~ .x*10.16)),
         across(c(Estimate,  Q2.5,   Q97.5), ~ round(.x, 2))) %>% 
  mutate(res = paste0(Estimate, " (", Q2.5, " to ", Q97.5, ")")) %>% 
  select(outcome, model_type, adj, res)

## extract priors ----
priors_smry <- map(mdls, ~ prior_summary(.x) %>% as_tibble()) %>% 
  bind_rows(.id = "out_term_type") %>% 
  SepModelNames() %>% 
  filter(!prior == "") %>% 
  select(-(coef:source)) %>% 
  spread(class, prior)
## Extract P values ----
ps <- map(mdls, function(x) {
  x <- as_draws_df(x, variable = "b_Intercept")$b_Intercept
  # print(any(p>=0))
  tibble(p = mean(x < 0), n = length(x))
}) %>% 
  bind_rows(.id = "out_term_type") %>% 
  SepModelNames() %>% 
  filter(term %in% c("como", "nter"))
ps <- ps %>% 
  mutate(p = if_else(p == 1, ">0.999", as.character(round(p, 3))))
t2 <- t2 %>% 
  inner_join(ps %>% filter(term == "como")) 
t2b <- t2 %>% 
  mutate(res = paste0(res, "; P = ", p)) %>% 
  select(outcome, model_type, adj, res) %>% 
  distinct() %>% 
  spread(outcome, res)

write_csv(ps, "Outputs/bayesianPs.csv")
write_csv(t2b, "Outputs/t2.csv")
write_csv(t2_unscale, "Outputs/t2_unscaled.csv")
write_csv(t3, "Outputs/t3.csv")



## extract condition and drug level effects ----
LinCombCond <- function(modelobj, df){
  recover_types(modelobj, df) %>%
    spread_draws(b_Intercept, r_sail[sail,]) %>%
    mean_qi(m = b_Intercept + r_sail, .width = c(0.95, 0.8, 0.5))
}
LinCombDrug <- function(modelobj, df){
 df <- recover_types(modelobj, df)
df %>% spread_draws(b_Intercept, r_code5[code5,]) %>%
       mean_qi(m = b_Intercept + r_code5, .width = c(0.95, 0.8, 0.5))
}

## loop through objects
modelnames <- c("eq5d_como_drugcond_adj", "mcs_como_drugcond_adj", "pcs_como_drugcond_adj", 
                "eq5d_nter_drugcond_adj", "mcs_nter_drugcond_adj", "pcs_nter_drugcond_adj")
dfnames   <- c("df_eq5d_como", "df_mcs_como", "df_pcs_como",
               "df_eq5d_nter", "df_mcs_nter", "df_pcs_nter")
cond_lvl <- map2(mdls[modelnames], all_df[dfnames], ~ LinCombCond(.x, .y)) %>% 
  bind_rows(.id = "out_term_type")  %>% 
  SepModelNames() %>% 
  mutate(term_l = factor(term, levels = c("como", "nter"), c("Comorbidity", "Trt-comorbidity interaction")),
         outcome_l = factor(outcome, levels = c("eq5d", "mcs", "pcs"), labels = c("EQ5D", "SF-36 PCS", "SF-36 MCS")),
         sail=replace(sail, sail=="Dementia (any)", "Dementia"), 
         sail=replace(sail, sail=="Parkinson's disease (all)", "Parkinson's Disease"))
drug_lvl <- map2(mdls[modelnames], all_df[dfnames], ~ LinCombDrug(.x, .y)) %>% 
  bind_rows(.id = "out_term_type")  %>% 
  SepModelNames() %>% 
  mutate(term_l = factor(term, levels = c("como", "nter"), c("Comorbidity", "Trt-comorbidity interaction")),
         outcome_l = factor(outcome, levels = c("eq5d", "mcs", "pcs"), labels = c("EQ5D", "SF-36 PCS", "SF-36 MCS"))) 

drug_lvl_for_plot <- drug_lvl %>% 
  mutate(code5 = replace(code5, code5 == "A10BG", "A10BG- Thiazolidinediones"), 
         code5 = replace(code5, code5 == "A10BJ", "A10BJ- Glucagon-Like Peptide-1 Analogues"),
         code5 = replace(code5, code5 == "A10BK", "A10BK- Sodium Glucose Co-transporter 2 Inhibitors"),
         code5 = replace(code5, code5 == "B01AE", "B01AE- Direct Thrombin Inhibitors"),
         code5 = replace(code5, code5 == "G04BE", "G04BE- Drugs Used in Erectile Dysfunction"),
         code5 = replace(code5, code5 == "H05AA", "H05AA- Parathyroid Hormones and Analogues"),
         code5 = replace(code5, code5 == "L01XE", "L01XE- Protein Kinase Inhibitors"),
         code5 = replace(code5, code5 == "L04AA", "L04AA- Selective Immunosuppressants"),
         code5 = replace(code5, code5 == "L04AB", "L04AB- Tumour Necrosis Factor Alpha Inhibitors"),
         code5 = replace(code5, code5 == "L04AC", "L04AC- Interleukin Inhibitors"),
         code5 = replace(code5, code5 == "M05BA", "M05BA- Bisphosphonates"),
         code5 = replace(code5, code5 == "N03AX", "N03AX- Other Antiepileptics"),
         code5 = replace(code5, code5 == "N04BC", "N04BC- Dopamine Agonists"),
         code5 = replace(code5, code5 == "R03DX", "R03DX- Other Systemic Drugs For Obstructive Airway Diseases"))


# This will draw a plot with distributions by condition
plot1 <- ggplot(cond_lvl,
                aes(y=reorder(sail, desc(sail)), x = m, xmin = .lower, xmax = .upper, colour = term_l)) +
  geom_pointinterval(position = position_dodge(0.5)) + 
  geom_vline(xintercept = c(0),
             linetype = c("dotted"), color = c("black")) +
  theme_classic() + 
  labs(x = "Effect estimate (50%, 80% and 95% credibility intervals)", y = "Condition") +
  facet_wrap( ~ outcome_l, scales = "free_x") +
  scale_color_colorblind("")+
  theme(text = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(angle = 30, hjust = 1), 
        panel.spacing.x = unit(0.7, "lines") ) 
plot1
tiff("Outputs/figure1.tiff", height = 1700*2, width = 2250*2, res = 600, units = "px", compression = "lzw")
plot1
dev.off()

# This will draw a plot with distributions by class
plot2 <- ggplot(drug_lvl_for_plot,
                aes(y=reorder(code5, desc(code5)), x = m, xmin = .lower, xmax = .upper, colour = term_l)) +
  geom_pointinterval(position = position_dodge(0.5)) + 
  geom_vline(xintercept = c(0),
             linetype = c("dotted"), color = c("black")) +
  theme_classic() + 
  labs(x = "Effect estimate (50%, 80% and 95% credibility intervals)", y = "Drug Treatment Comparisons") +
  facet_wrap(~ outcome_l, scales = "free_x") +
  scale_color_colorblind("")+
  theme(text = element_text(size = 8, family = "serif"),
        axis.text.x = element_text(angle =30, hjust = 1))

plot2
tiff("Outputs/figure2.tiff", height = 1700*2, width = 2250*2, 
     res = 600, units = "px", compression = "lzw")
plot2
dev.off()
## predict posterior; create priors for subsequent analyses
## All index conditions for
## 1. drug-classes we observe for that index condition
## 2. Some unobserved drug class we have never seen
## 3. Some unobserved condition we have never seen

## Most convenient way to create the estimates for the unknown condition and class is to 
## set the standard error for an individual trial to zero, gives same result as sampling from posterior and calculating
## for each dataset create lookup variables for the sail variable and drug class variable
Createlkp <- function(mydf = all_df$df_eq5d_nter, mydf_comp =all_df_cmp$df_eq5d_nter){
  a <- bind_cols(mydf %>% count(code5) %>% rename(lbl_cnt = n, lbl = code5), (tibble(lvl = mydf_comp$code5) %>% count(lvl)) %>% rename(lvl_cnt = n))
  b <- bind_cols(mydf %>% count(sail)  %>% rename(lbl_cnt = n, lbl = sail),  (tibble(lvl = mydf_comp$sail)  %>% count(lvl)) %>% rename(lvl_cnt = n))
  ab <- bind_rows(code5 = a,
                  sail = b)
  if(!all(ab$lbl_cnt == ab$lvl_cnt)) warning ("label mismatches")
  sail <- b$lvl
  names(sail) <- b$lbl
  code5 <- a$lvl
  names(code5) <- a$lbl
  ab
  list(sail = sail, code5 = code5)
}
eq5d_lkp <- Createlkp()
mcs_lkp  <- Createlkp(all_df$df_mcs_nter, all_df_cmp$df_mcs_nter)
pcs_lkp  <- Createlkp(all_df$df_pcs_nter, all_df_cmp$df_pcs_nter)

## Produce summary stats for each prior AND diagnostic plot of the fit of this
ObtainPosteriors <- function(mydf, mdl_choose = mdls$eq5d_nter_drugcond, mylkp){
  choose_outcome <- mydf %>% 
    distinct(sail, code5) %>% 
    rename(sail_lbl = sail, code5_lbl = code5) %>% 
    mutate(sail = mylkp$sail[sail_lbl],
           code5 = mylkp$code5[code5_lbl],
           std.error = 0) %>% 
    arrange(sail_lbl, code5_lbl)

  ## Create new unknown drug class for each condition
  choose_outcome2 <- choose_outcome %>%
    distinct(sail, .keep_all = TRUE) %>%
    mutate(code5_lbl = "Unknown", code5 = max(code5 + 100))
  ## Create new unknown drug class and unknown condition condition
  choose_outcome3 <- choose_outcome2 %>% 
    distinct(code5_lbl, .keep_all = T) %>% 
    mutate(sail_lbl = "Unknown", sail = max(choose_outcome$sail) + 900)
  
  ## Join to original observed data
  choose_outcome <- bind_rows(choose_outcome,
                         choose_outcome2,
                         choose_outcome3)
  
  choose_outcome$res <- posterior_predict(mdl_choose,
                    newdata = choose_outcome,
                    re_formula = ~ (1|sail) + (1|code5),
                    allow_new_levels = TRUE,
                    # sample_new_levels = "gaussian"
                    ) %>% 
    as.data.frame() %>% 
    as.list()
  choose_outcome$res_m <- map_dbl(choose_outcome$res, mean)
  choose_outcome$res_s <- map_dbl(choose_outcome$res, sd)
  
  ## examine distributions
  choose_outcome$res <- map(choose_outcome$res, sample, size = 1000)
  choose_outcome$res <- map(choose_outcome$res, ~ tibble(y = .x))
  choose_outcome
}
all_post <- bind_rows(
  eq5d = ObtainPosteriors(all_df$df_eq5d_nter,mdl_choose = mdls$eq5d_nter_drugcond, mylkp = eq5d_lkp),
  mcs  = ObtainPosteriors(all_df$df_mcs_nter, mdl_choose = mdls$mcs_nter_drugcond, mylkp = mcs_lkp),
  pcs  = ObtainPosteriors(all_df$df_pcs_nter, mdl_choose = mdls$pcs_nter_drugcond, mylkp = pcs_lkp),
                          .id = "outcome")

## Summarise posteriors using t-distribution; this bit is slow as running for all combinations ----
mdl_create <- brm(y ~ 1, data = all_post$res[[1]], family = "student", chains = 0)
all_post$t_mdl <- map(all_post$res, ~ update(mdl_create, recompile = FALSE, chains = 4, newdata = .x))
all_post$t_smry <- map(all_post$t_mdl, function(mdl) {
    a <- summary(mdl)
    tibble(m = a$fixed[1,1], s = a$spec_pars[1,1], df = a$spec_pars[2,1])
  })

## Sample from t-distribution to compare to raw samples and plot
all_post$t_smpl <- map(all_post$t_smry, ~ rt(1000, df=.x$df)*.x$s + .x$m)
all_post_lng <- all_post %>% 
  select(outcome, sail_lbl, code5_lbl, res, t_smpl) %>% 
  unnest(cols = c(res, t_smpl)) %>% 
  rename(r_smpl = y) %>% 
  gather("smpl", "value", r_smpl, t_smpl)

plot1 <- ggplot(all_post_lng, aes(x = value, colour = smpl)) + 
  geom_density(alpha = 0.2, fill = NA, size = 0.1) + 
  facet_wrap(outcome ~ sail_lbl + code5_lbl)
pdf("Outputs/summarise_prior_t_distribution.pdf", height = 20, width = 40)
plot1
dev.off()

svg("Outputs/summarise_prior_t_distribution.svg", height = 20, width = 40)
plot1
dev.off()

## Summarise priors for appendix
all_post_prnt <- all_post %>% 
  select(outcome, sail_lbl, code5_lbl, t_smry) %>% 
  unnest(t_smry) 
write_csv(all_post_prnt, "Outputs/summarise_posteriors.csv")
