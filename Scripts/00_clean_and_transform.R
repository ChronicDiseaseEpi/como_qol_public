# 00_clean_and_transform.R
library(tidyverse)
library(brms)
library(tidybayes)

## Functions ----
SepModelNames <- function(x) { x %>% 
    separate(out_term_type, c("outcome", "term", "model_type", "adj"), sep = "_")
  }

## Load data ----
cfs <- read_csv("Data/cfs.csv")

## List of models ----
cept      <- c(
  # age and sex adjusted, no impler models
  m5 = "(Intercept), base, age_std, sex; (Intercept), base, age_std")
arm       <- c(
  # no arm interactions, just age-sex adjusted
  m54 = "(Intercept), base, sex, arm, arm:sex; (Intercept), base, arm")
como_count <- c(
  # Unadjusted
  m21 = "(Intercept), base, como_cnt",
  # Age adjusted
  m25 = "(Intercept), base, age_std, como_cnt", 
  # Age and sex adjusted
  m45 = "(Intercept), base, sex, como_cnt; (Intercept), base, como_cnt",
  # age and sex adjusted
  m1  = "(Intercept), base, age_std, sex, como_cnt")

como_nter <- c(
  # no other interactions
  m23 = "(Intercept), base, como_cnt, arm, arm:como_cnt",
  # also age interaction
  m32 = "(Intercept), base, arm, age_std, como_cnt, arm:age_std, arm:como_cnt", 
  # also sex interaction
  m52 = "(Intercept), base, sex, arm, como_cnt, arm:sex, arm:como_cnt; (Intercept), base, arm, como_cnt, arm:como_cnt",
  # also age and sex interaction
  m8 = "(Intercept), base, arm, age_std, sex, como_cnt, arm:age_std, arm:sex, arm:como_cnt;(Intercept), base, arm, age_std, como_cnt, arm:age_std, arm:como_cnt; (Intercept), base, arm, age_std, sex, como_cnt, arm:age_std, arm:sex")

## Identify outcomes ----
myoutcome  <- "eq5d"
myoutcome2 <- "sf36mcs"
myoutcome3 <- "sf36pcs"

## identify is positive trial using simplest model with main arn effects
m54 <- cfs %>% 
  filter(model_n == "m54", term == "arm") %>% 
  mutate(pstv = if_else( (estimate - 1.96*std.error) > 0, 1L, 0L))
postrial <- m54 %>% 
  filter(pstv == 1) %>% 
  distinct(outcome, nct_id)
postrial %>% 
  count(outcome)
negtrial <-  m54 %>% 
  filter(pstv == 0) %>% 
  distinct(outcome, nct_id)
negtrial %>% 
  count(outcome)
rm(m54, negtrial)

## compose data using tidybayes ----
## want to divide estimate and SE by SD so can use same priors on all models
all_df <- list(
  df_eq5d_unad = identity(cfs %>% 
                            filter(model_n == "m21",
                                   outcome == "eq5d",
                                   term == "como_cnt")),
  df_eq5d_nter = identity(cfs %>% 
                               filter(model_n == "m8",
                                      outcome == "eq5d",
                                      term == "arm:como_cnt")),
  df_eq5d_como = identity(cfs %>% 
                                 filter(model_n == "m1",
                                        outcome == "eq5d",
                                        term == "como_cnt")),
  df_mcs_unad = identity(cfs %>% 
                           filter(model_n == "m21",
                                  outcome == "sf36mcs",
                                  term == "como_cnt")),
  
  df_mcs_nter = identity(cfs %>% 
                                filter(model_n == "m8",
                                       outcome == "sf36mcs",
                                       term == "arm:como_cnt")),
  
  df_mcs_como = identity(cfs %>% 
                                filter(model_n == "m1",
                                       outcome == "sf36mcs",
                                       term == "como_cnt")),
  df_pcs_unad = identity(cfs %>% 
                           filter(model_n == "m21",
                                  outcome == "sf36pcs",
                                  term == "como_cnt")),
  df_pcs_nter = identity(cfs %>% 
                                filter(model_n == "m8",
                                       outcome == "sf36pcs",
                                       term == "arm:como_cnt")),
  
  df_pcs_como = identity(cfs %>% 
                                filter(model_n == "m1",
                                       outcome == "sf36pcs",
                                       term == "como_cnt")),
  
  df_eq5d_comosq = identity(cfs %>% 
                                   filter(model_n == "m2",
                                          outcome == "eq5d",
                                          term == "I(como_cnt^2)")),
  df_mcs_comosq = identity(cfs %>% 
                                  filter(model_n == "m2",
                                         outcome == "sf36mcs",
                                         term == "I(como_cnt^2)")),
  df_pcs_comosq = identity(cfs %>% 
                                  filter(model_n == "m2",
                                         outcome == "sf36pcs",
                                         term == "I(como_cnt^2)"))

)
all_df <- map(all_df, ~ .x %>% 
                select(estimate, std.error, outcome, sail, code5, nct_id))
all_df[c("df_eq5d_unad", "df_eq5d_nter", "df_eq5d_como", "df_eq5d_comosq")] <-
  map(all_df[c("df_eq5d_unad","df_eq5d_nter", "df_eq5d_como", "df_eq5d_comosq")], ~ .x %>% 
          mutate_at(vars(estimate, std.error), function(x) x / 0.23))
all_df[c("df_mcs_unad","df_mcs_nter", "df_mcs_como", "df_mcs_comosq")] <- 
  map(all_df[c("df_mcs_unad", "df_mcs_nter", "df_mcs_como", "df_mcs_comosq")], ~ .x %>% 
          mutate_at(vars(estimate, std.error), function(x) x / 9.08))
all_df[c("df_pcs_unad","df_pcs_nter", "df_pcs_como", "df_pcs_comosq")] <-
  map(all_df[c("df_pcs_unad","df_pcs_nter", "df_pcs_como", "df_pcs_comosq")], ~ .x %>% 
          mutate_at(vars(estimate, std.error), function(x) x / 10.16))
## Convert data into lists for passing to model
all_df_cmp <- map(all_df, tidybayes::compose_data)
all_df_cmp_pstv <- map(all_df, ~ tidybayes::compose_data(.x %>% semi_join(postrial)))

## Draw plot of raw effect estimates
all_df_lng <- bind_rows(all_df, .id = "df")
all_df_lng <- all_df_lng %>% 
  arrange(sail, code5) 
lvl_assign <- unique(all_df_lng$nct_id)
all_df_lng <- all_df_lng %>% 
  mutate(nct_id = factor(nct_id, levels = lvl_assign, ordered = TRUE))

plot1 <- ggplot(all_df_lng, aes(x = nct_id, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, colour = interaction(sail, code5))) +
  geom_point() + 
  geom_linerange() +
  facet_wrap(~ df, scales = "free") +
  coord_flip()
pdf("Outputs/raw_effect_estimates.pdf", width = 20, height = 15)
plot1
dev.off()
rm(all_df_lng, plot1)
