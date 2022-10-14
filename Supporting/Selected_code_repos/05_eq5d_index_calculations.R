##Further harmonisation of EQ5D and calculation of index scores
library(tidyverse)
library(readr)

#Read in required dfs

eq5d <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Processed_data/eq5d.Rds")

new_labels <- read_delim("E:/Research Project 1732/Consolidate_outcomes/Created_metadata/resolve_eq5d_question_names.csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)

eq5d5_index <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Created_metadata/eq5d5.Rds")
eq5d3_index <- readRDS("E:/Research Project 1732/Consolidate_outcomes/Created_metadata/eq5d3.Rds")

#harmonise parameter labels
eq5d_labelled <- eq5d %>% 
  inner_join(new_labels) %>% 
  select(condition, nct_id, subjid, visit, parameter = assigned, value)

##Check to see if any not labelled

x <- eq5d %>% anti_join(new_labels)
table(x$question)

y <- eq5d %>% filter(question == "")
table(y$value)
   #therefore shows that any rows where question is.na, value is.na therefore ok to be missing from eq5d_labelled

#Check which version of Eq5d for each trial (3L vs 5L)- i.e. is the max score 3 or 5

#remove visual analog score and provided index scores for this exercise

#These are the three trials which have values 4 or 5 therefore must be from eq5d-5l

#Remaining 31 trials- documentation reviewed to ensure EQ5D-3L- all EQ5D_3l

#Join with index scores to produce final 3l and 5l dataframes
eq5d_5l <- eq5d_labelled %>% filter(nct_id == "NCT01474512" | nct_id == "NCT01597245" | nct_id == "NCT01646177")
eq5d_3l <- eq5d_labelled %>%
  filter(!nct_id == "NCT01474512") %>%
  filter(!nct_id == "NCT01597245") %>%
  filter(!nct_id == "NCT01646177")

eq5d_5l_final <- eq5d_5l %>% 
  group_by_at(vars(-value)) %>% 
  mutate(row_id = 1:n()) %>% ungroup() %>% 
  spread("parameter", "value") %>% 
  select(condition, nct_id, subjid, visit, MO, SC, UA, PD, AD, Provided_index_UK, VAS) %>% 
  mutate_at(c(5,6,7,8,9,10,11), as.numeric) %>% 
  inner_join(eq5d5_index) %>% 
  select(condition, nct_id, subjid, visit, MO, SC, UA, PD, AD, calculated_index = score, Provided_index_UK, VAS )
# write_rds(eq5d_5l_final, "E:/Research Project 1732/Consolidate_outcomes/Processed_data/eq5d_5l_final")

eq5d_3l_final <- eq5d_3l %>% 
  group_by_at(vars(-value)) %>% 
  mutate(row_id = 1:n()) %>% ungroup() %>% 
  spread("parameter", "value") %>% 
  select(condition, nct_id, subjid, visit, MO, SC, UA, PD, AD, Provided_index_UK, VAS) %>% 
  mutate_at(c(5,6,7,8,9,10,11), as.numeric) %>% 
  inner_join(eq5d3_index) %>% 
  select(condition, nct_id, subjid, visit, MO, SC, UA, PD, AD, calculated_index = score, Provided_index_UK, VAS )

#Join 3l and 5l together as index now  comparable
eq5d_final <- bind_rows(eq5d_5l_final,
                        eq5d_3l_final)
## Note two missing trials 

write_rds(eq5d_final, "E:/Research Project 1732/Consolidate_outcomes/Processed_data/eq5d_score_calculated.Rds")
