#05_plot_trial_diffs
# library(tidyverse)
source("Scripts/01a_baseline_como_qol.R")

## treatment arm only
m54 <- cfs %>% 
  filter(model_n == "m54",
         term %in% c("arm", "arm:sex", "sex", "(Intercept)", "base"))
m54 <- m54 %>% 
  select(-std.error) %>% 
  spread(term, estimate, fill = 0)

base_smry_fig2 <- base_smry_fig %>% 
  rename(outcome = parameter) %>% 
  mutate(outcome = if_else(outcome == "eq5d_index", "eq5d", outcome) %>% str_remove("_")) %>% 
  select(-outcome_s) %>% 
  mutate(arm = if_else(arm == 1, "trtd", "cntl"),
         sex = if_else(sex ==1, "men_pre", "wom_pre")) %>% 
  pivot_wider(names_from = c(arm, sex), values_from = outcome_m)

m54 <- m54 %>% 
  left_join(base_smry_fig2)

m54 <- m54 %>% 
  ## drop this step as yoda data exported
  # group_by(outcome) %>%
  # mutate_at(vars(cntl_wom_pre, cntl_men_pre, trtd_wom_pre, trtd_men_pre), function(x) if_else(is.na(x), mean(x, na.rm = TRUE), x)) %>%
  ungroup() %>%
  ## 
  mutate(cntl_wom_pos = `(Intercept)` + (base-1)*cntl_wom_pre       + cntl_wom_pre,
         cntl_men_pos = `(Intercept)` + (base-1)*cntl_men_pre + sex + cntl_men_pre,
         trtd_wom_pos = `(Intercept)` + (base-1)*trtd_wom_pre + arm + trtd_wom_pre,
         trtd_men_pos = `(Intercept)` + (base-1)*trtd_men_pre + 
           sex + arm + `arm:sex`                                  + trtd_men_pre) %>% 
  select(-arm, -`arm:sex`, -sex, -`(Intercept)`, -base)

## Note, expecting to drop men and women from "NCT00670501" and "NCT00439647" respectively
m54 <- m54 %>%
  gather("alloc_sex", "estimate", 
         cntl_wom_pre, cntl_men_pre, trtd_wom_pre, trtd_men_pre,
         cntl_wom_pos, cntl_men_pos, trtd_wom_pos, trtd_men_pos, na.rm = TRUE)
m54 <- m54 %>%
  separate(alloc_sex, into = c("alloc", "sex", "time"), sep = "_") 

## Summary statistics for whole study level ----
m54_smry <- m54 %>% 
  select(repo, nct_id, outcome, alloc, sex, time, estimate) %>% 
  spread(time, estimate)
base_smry_text2 <- base_smry_text %>% 
  mutate(parameter = parameter %>% 
           str_remove("index") %>% 
           str_remove("_"))
m54_smry <-  m54_smry %>% 
               mutate(sex = if_else(sex == "men", 1L, 0L),
                      arm = if_else(alloc == "cntl", 0L, 1L)) %>% 
  select(repo, nct_id, parameter = outcome, arm, sex, pos, pre )
m54_smry <- m54_smry %>% 
  inner_join(base_smry_text2)
m54_smry_arm <- m54_smry %>% 
  group_by(parameter, arm) %>% 
  summarise(across(c(pre, pos), ~ weighted.mean(.x, n))) %>% 
  ungroup() %>% 
  mutate(across(c(pre, pos), ~ if_else(parameter == "eq5d", round(100*.x), round(.x))),
         parameter = if_else(parameter == "eq5d", "eq5dx100", parameter))
m54_smry_arm_yoda <- m54_smry %>% 
  filter(repo  == "yoda") %>% 
  group_by(nct_id, parameter) %>% 
  summarise(across(c(pre, pos), ~ weighted.mean(.x, n))) %>% 
  ungroup()

write_csv(m54_smry_arm, "Outputs/raw_change_qol.csv")

## Group for plotting ----
m54 <- m54 %>% 
  group_by(nct_id, outcome, sex, alloc) %>% 
  mutate(imprv = if_else(estimate[time == "pos"] > 
                           estimate[time == "pre"] , "black", "blue"),
         pre_tx_pl = case_when(
           time == "pre" ~ "black",
           time == "pos" & alloc == "cntl" ~ "blue",
           time == "pos" & alloc == "trtd" ~ "red")) %>% 
  ungroup() %>% 
  mutate(sex_l = if_else(sex == "men", "Men", "Women"),
         outcome_l = factor(outcome, levels = c("eq5d", "sf36mcs", "sf36pcs"),
                            labels = c("EQ5D", "SF-36 MCS", "SF-36 PCS")),
         alloc_l = factor(alloc,
                          levels = c("cntl",
                                     "trtd"),
                          labels = c("Comparator",
                                     "Treatment")),
         time_l = factor(time, levels = c("pre", "pos"), c("Baseline", "Follow-up")),
         xvar   = paste(alloc_l, time_l, sep = "\n"))

## Apendix plot for change pre and pot in treatmetn and control arms
plot2 <- ggplot(m54, aes(x = xvar, y = estimate, colour = imprv, group = interaction(nct_id, alloc_l))) +
  geom_point() +
  geom_line(alpha = 0.5, mapping = aes(colour = imprv)) +
  facet_wrap(outcome_l ~ sex_l, scales = "free_y", ncol = 2) +
  scale_colour_identity() +
  scale_y_continuous("Quality of life score") +
  scale_x_discrete("") +
  theme(text = element_text(size = 8))
plot2

m54 %>% 
  group_by(outcome, sex, alloc_l, time_l) %>% 
  summarise(m = mean(estimate)) %>% 
  ungroup() %>% 
  spread(time_l, m) %>% 
  mutate(diff = `Follow-up` - Baseline)

## compare magnitude of treatment effect to magnitude of comorbidity effect
# m8
m3 <- cfs %>% 
  filter(model_n == "m3",
         term %in% c("como_cnt", "arm")) %>% 
  select(-std.error) %>% 
  pivot_wider(names_from = term, values_from = c(estimate)) %>% 
  mutate(outcome_l = factor(outcome, levels = c("eq5d", "sf36mcs", "sf36pcs"),
                            labels = c("EQ5D", "SF-36 MCS", "SF-36 PCS")))

plot3 <- ggplot(m3, aes(x = arm, y = como_cnt)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~outcome_l, scales = "free") +
  scale_y_continuous("Comorbidity count association") +
  scale_x_continuous("Treatment effect") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Comparison of magnitude of effect of treatment versus magnitude of comorbidity association")
plot3

plot3b <- ggplot(m3 %>% filter(!nct_id == "NCT01191268"), aes(x = outcome_l, y = como_cnt/arm)) +
  geom_violin() +
  # facet_wrap(~outcome_l, scales = "free") +
  # scale_x_continuous("Comorbidity count/treatment effect ratio") +
  ggtitle("Comparison of magnitude of effect of treatment versus magnitude of comorbidity association")
plot3b


## compare magnitude of treatment effect to magnitude of treatment-comorbidity interaction
# m8
m3 <- cfs %>% 
  filter(model_n == "m3",
         term %in% c("arm", "arm:como_cnt")) %>% 
  select(-std.error) %>% 
  pivot_wider(names_from = term, values_from = c(estimate)) %>% 
  mutate(outcome_l = factor(outcome, levels = c("eq5d", "sf36mcs", "sf36pcs"),
                            labels = c("EQ5D", "SF-36 MCS", "SF-36 PCS")))
plot4 <- ggplot(m3, aes(x = arm, y = `arm:como_cnt`)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~outcome_l, scales = "free") +
  # scale_x_continuous("Comorbidity count association") +
  # scale_y_continuous("Treatment effect") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Comparison of magnitude of effect of treatment versus magnitude of comorbidity treatment interaction")
plot4

plot4b <- ggplot(m3 %>% filter(!nct_id == "NCT01191268"), aes(x = outcome_l, y = `arm:como_cnt`/arm)) +
  geom_violin() +
  # facet_wrap(~outcome_l, scales = "free") +
  # scale_x_continuous("Comorbidity count/treatment effect ratio") +
  ggtitle("Comparison of magnitude of effect of treatment versus magnitude of comorbidity treatment interaction")
plot4b
