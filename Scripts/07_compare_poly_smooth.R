# compare_poly_smooth
library(tidyverse)
## Note that some standard errors are zero due to rounding. They are not truly zero
mydf <- read_csv("Data/qol_models_compare_sq_sm.csv")
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
  facet_wrap(~nct_id + outcome, scales = "free_y") +
  scale_y_continuous("") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_colour_discrete("")
plot1
png("Outputs/compare_poly_sq.png", height = 12, width = 16, unit = "in", res = 300)
plot1
dev.off()
