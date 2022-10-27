# compare_poly_smooth
library(tidyverse)
## Note that some standard errors are zero due to rounding. They are not truly zero
mydf <- read_csv("Data/qol_models_compare_sq_sm.csv")
mydf2 <- read_csv("Data/qol_models_compare_sq_sm_csdr.csv")
mydf2 <- mydf2 %>% 
  rename(modeltype = model_type) %>% 
  mutate(modeltype = if_else(modeltype == "sp", "comosm", "comosp")) %>% 
  filter(term == "como") %>% 
  select(-term)
mydf <- mydf %>% 
  bind_rows(mydf2)
rm(mydf2)

mydf <- mydf %>% 
  gather("como_cnt", "res", `0`:`10`, na.rm = TRUE) %>% 
  separate(res, into = c("est", "se"),sep = "\\(") %>% 
  mutate(se = str_remove(se, "\\)") %>% as.double(),
         est = as.double(est),
         como_cnt = as.integer(como_cnt),
         lci = est - 1.96*se,
         uci = est + 1.96*se,
         modeltype = if_else(modeltype == "comosm", "Pen. splines", "Polynomial"))
plot1 <- ggplot(mydf, aes(x = como_cnt, y = est, ymin = lci, ymax = uci, colour = modeltype)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_wrap(~nct_id, scales = "free_y") +
  scale_y_continuous("") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_colour_discrete("")
plot1

png("Outputs/compare_poly_sq_eq5d.png", height = 12, width = 16, unit = "in", res = 300)
plot1 %+% (mydf %>% filter(outcome == "eq5d_index"))
dev.off()

png("Outputs/compare_poly_sq_pcs.png", height = 12, width = 16, unit = "in", res = 300)
plot1 %+% (mydf %>% filter(outcome %in% c("sf36_pcs", "sf36_PCS")))
dev.off()

png("Outputs/compare_poly_sq_mcs.png", height = 12, width = 16, unit = "in", res = 300)
plot1 %+% (mydf %>% filter(outcome %in% c("sf36_mcs", "sf36_MCS")))
dev.off()
