##Consolidate SF36 scores

library(tidyverse)
library(haven)
library(readr)
library(readxl)
library(tidyr)

#Code and explanation from github saved in ".../Consolidate_outcomes/Scratch_data/SF-36 Scoring.Rmd"

#Read in SF36 table #----

sf36 <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Processed_data/sf36.Rds")


#Standardise domain labels #----
#Not all trials have total or components scores, remove for now. 
#filter to remove entries with blank visit name as will not be able to distinguish
sf36 <- sf36 %>%
  mutate(question = replace(question, question == "bp" | question == "Bodily Pain" | question == "SF36 Bodily pain", "BP")) %>%
  mutate(question = replace(question, question == "gh" | question == "General Health" | question == "SF36 General health", "GH")) %>%
  mutate(question = replace(question, question == "mh" | question == "Mental Health" | question == "SF36 Mental health", "MH")) %>%
  mutate(question = replace(question, question == "pf" | question == "Physical Functioning" | question == "SF36 Physical functioning", "PF")) %>%
  mutate(question = replace(question, question == "re" | question == "Role Emotional" | question == "SF36 Role limitations emotional problems", "RE")) %>%
  mutate(question = replace(question, question == "rp" | question == "Role Physical" | question == "SF36 Role limitations physical problems", "RP")) %>%
  mutate(question = replace(question, question == "sf" | question == "Social Functioning" | question == "SF36 Social functioning", "SF")) %>%
  mutate(question = replace(question, question == "vt" | question == "Vitality" | question == "SF36 Vitality", "VT")) %>%
  filter(!question == "total") %>%
  filter(!question == "SF36 Physical component summary") %>%
  filter(!question == "SF36 Mental component summary") %>%
  mutate(value = as.numeric(value)) %>%
  filter(!visit == "") %>%
  distinct()



#On review the following 2 trials are already standardised z scores #----
#NCT00783718 and NCT00224171 therefore do not standardise again

#Bring in General US population values #----

SF36_GenPopUS <- read_excel("Scratch_data/Copy of SF36_GenPopUS_population_median_sd.xlsx")

#Spread domain scores per person per visit for code to work #----
subScale8 <- sf36 %>% 
  group_by(condition, nct_id, subjid, question, visit) %>% 
  summarise(value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = question, values_from = value)

## calculate standardized subscale scores for each trial participant per visit #----
std_subScale8 <- subScale8 %>% 
  # filter(!nct_id== "NCT00783718") %>% 
  # filter(!nct_id== "NCT01224171") %>% 
  mutate(
  
    # divide sample subscale score by the mean scale score from general US population and divide by std dev from general US population
    stdPF = (PF - as.numeric(SF36_GenPopUS[1,3])) / as.numeric(SF36_GenPopUS[1,5]),
    stdRP = (RP - as.numeric(SF36_GenPopUS[2,3])) / as.numeric(SF36_GenPopUS[2,5]),
    stdBP = (BP - as.numeric(SF36_GenPopUS[3,3])) / as.numeric(SF36_GenPopUS[3,5]),
    stdGH = (GH - as.numeric(SF36_GenPopUS[4,3])) / as.numeric(SF36_GenPopUS[4,5]),
    stdVT = (VT - as.numeric(SF36_GenPopUS[5,3])) / as.numeric(SF36_GenPopUS[5,5]),
    stdSF = (SF - as.numeric(SF36_GenPopUS[6,3])) / as.numeric(SF36_GenPopUS[6,5]),
    stdRE = (RE - as.numeric(SF36_GenPopUS[7,3])) / as.numeric(SF36_GenPopUS[7,5]),
    stdMH = (MH - as.numeric(SF36_GenPopUS[8,3])) / as.numeric(SF36_GenPopUS[8,5])
   )


std_subScale8_lng <- std_subScale8 %>% 
  select(nct_id, subjid, visit, stdPF, stdRP, stdBP, stdGH, stdVT, stdSF, stdRE, stdMH) %>% 
  gather("question", "value",   stdPF, stdRP, stdBP, stdGH, stdVT, stdSF, stdRE, stdMH, na.rm = TRUE)

multiplier <- tibble(question = c("stdPF", "stdRP", "stdBP", "stdGH", "stdVT", "stdSF", "stdRE", "stdMH"),
                     wting  = c(10, 4, 2, 5, 4, 2, 3, 5),
                     type = c("p", "p", "p", "p", "m", "m", "m", "m"))

std_subScale8_lng <- std_subScale8_lng %>% 
  inner_join(multiplier) %>% 
  mutate(value_mult = value * wting) %>% 
  group_by(nct_id, subjid, visit, type) %>% 
  summarise(numer = sum(value_mult),
            denom = sum(wting)) %>% 
  ungroup() %>% 
  mutate(std_score = (numer/denom) * 10 + 50) %>% 
  select(-numer, -denom)

std_subScale8_lng <- std_subScale8_lng %>% 
  mutate(parameter = if_else(type == "p", "sf36_pcs", "sf36_mcs")) %>% 
  select(-type) %>% 
  rename(value = std_score)

saveRDS(std_subScale8_lng, "Processed_data/sf36_score_calculated.Rds")

pdf("Outputs/check_sf_36_scores.pdf")
for (i in seq_along(std_subScale8_lng_nst$nct_id)) {
  plot_chk <- ggplot(std_subScale8_lng_nst$data[[i]], aes(x = std_score, colour = type, fill = type)) + 
    geom_histogram(alpha = 0.2) + 
    facet_wrap(~ visit) +
    ggtitle(std_subScale8_lng_nst$nct_id[i])
  invisible(print(plot_chk))
}
dev.off()



#Bind together EQ5d and SF36 final versions for modelling. #----

sf36_score_calculated <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Processed_data/sf36_score_calculated.Rds")

eq5d_score_calculated <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Processed_data/eq5d_score_calculated.Rds")

eq5d_score_calculated <- eq5d_score_calculated %>% 
  mutate(parameter = "eq5d_index") %>% 
  mutate(value = calculated_index) %>% 
  select(nct_id, subjid, visit, parameter, value)

#Bind together

eq5d_sf36_combined <- sf36_score_calculated %>% 
  bind_rows(eq5d_score_calculated)

saveRDS(eq5d_sf36_combined, "Processed_data/eq5d_sf36_combined.Rds")
