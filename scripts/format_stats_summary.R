# stats  model summary
library(dplyr)
library(stringr)
library(reshape2)

# read in csv with model results
setwd("C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc")
model_summary <- read.csv("output/formatted_stats_model_output_v3.csv")
# new formatting version:
# 0) select spatial extent?

# 1) select all rows that I want to transpose, keeping id rows (park, site, model, model#)
stats_estimates <- formatted_model_output %>%
  select(-delta_aicc, 
         -aicc, 
         -log_likelihood, 
         -adj_r2,
         -k)

# 2) transpose using cast function
stats_data_wide <- dcast(stats_estimates, extent+model+model_num~predictor, 
                          value.var = "estimate")

# 3) select delta aic, log likelihood and merge with the column above
stats_model_fit <- formatted_model_output %>% 
  select(extent, 
         model,
         model_num,
         k,
         delta_aicc, 
         log_likelihood, 
         adj_r2) %>%
  distinct()
  
merged_stats_data <- stats_model_fit %>%
  full_join(stats_data_wide, by=c("extent", "model", "model_num"))

write.csv(merged_stats_data, "output/formatted_stats_model_output_wide_dec14.csv", row.names=FALSE)
