# This script is used to compare general linear models of 7 key stream metrics 
# at three spatial extents using AICc (AIC for small sample sizes).

# Sections:
# 1) Read in and compile data for each spatial extent
# 2) Rescale the data from 0-1 so results can be compared between models
# 3) Create functions for each model comparison
# 4) Use functions to compare models at each spatial extent
# 5) Format and export model results 
# 6) Generate summary plot

# 1) Read in and compile data ----

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(kableExtra)
library(tibble)
library(stringr)
library(wiqid)
library(nlme)

# read in data ----
setwd("C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc")

disturbance_data_large <- read.csv("output/disturbance_data_large.csv")
disturbance_data_med <- read.csv("output/disturbance_data_med.csv")
disturbance_data_small <- read.csv("output/disturbance_data_small_proportional.csv")
field_data <- read.csv("output/all_site_data_means.csv")

# function to combine into one dataset and select relevant data
combine_data <- function(join_data){
  site_data <- field_data %>% 
    full_join(join_data, 
              by=c("park", "site")) %>% 
    select(park, 
           site, 
           catchment_area,
           total_disturbance,
           logging_percent, 
           insect_percent,
           total_road_density, 
           high_human_impact, 
           wetted_width, 
           wolmanD50, 
           percent_wetland, 
           percent_barren,
           depth, 
           flow, 
           percent_lake, 
           canopy,
           invertebrate_biomass, 
           tn, 
           conductivity, 
           embeddedness,
           ept_index, 
           periphyton_biomass, 
           shredder)
  
  return(site_data)
}

site_data_large <- combine_data(disturbance_data_large)
site_data_med <- combine_data(disturbance_data_med)
site_data_small <- combine_data(disturbance_data_small)

# 2) Rescale data with min max scaling ----

# Rescale function
min_max_scale <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

rescale_data <- function(data){
  rescaled_data <- data %>% 
    mutate(total_disturbance= min_max_scale(data$total_disturbance),
           total_road_density = min_max_scale(data$total_road_density),
           high_human_impact = min_max_scale(data$high_human_impact),
           wetted_width = min_max_scale(data$wetted_width),
           wolmanD50 = min_max_scale(data$wolmanD50),
           percent_wetland = min_max_scale(data$percent_wetland),
           depth = min_max_scale(data$depth),
           flow = min_max_scale(data$flow),
           percent_lake = min_max_scale(data$percent_lake),
           canopy = min_max_scale(data$canopy),
           invertebrate_biomass = min_max_scale(data$invertebrate_biomass),
           tn = min_max_scale(data$tn),
           conductivity = min_max_scale(data$conductivity),
           embeddedness = min_max_scale(data$embeddedness),
           ept_index = min_max_scale(data$ept_index),
           periphyton_biomass = min_max_scale(data$periphyton_biomass),
           percent_barren = min_max_scale(data$percent_barren),
           shredder = min_max_scale(data$shredder))
  
  return(rescaled_data)
}

rescaled_data_large <- rescale_data(site_data_large)
rescaled_data_med <- rescale_data(site_data_med)
rescaled_data_small <- rescale_data(site_data_small)

# 3) Functions for each model comparison ----

# function for extracting model results
extract_model_summary <- function(model, model_number, model_name){
  
  # calculate AIC
  aic <- wiqid::AICc(model)
  
  # extract items from the model output
  model_summary <- summary(model)
  
  r2 <- model_summary$adj.r.squared
  
  log_likelihood <- logLik(model)[1]
  
  # select coefficients and convert to dataframe
  mod_summary_coefficients <- model_summary$coefficients
  model_summary_df <- mod_summary_coefficients |> data.frame() |> 
    tibble::rownames_to_column(var = "predictor_var")
  
  # add a column for all values in summary table
  mod_summary <- model_summary_df |> 
    dplyr::mutate(model_num = model_number, 
                  model = model_name, 
                  adj_r2 = r2,
                  aic = aic,
                  log_like = log_likelihood) |> 
    dplyr::rename("estimate" = "Estimate", 
                  "std_error" = "Std..Error", 
                  "t_val" = "t.value", 
                  "p_val" = "Pr...t..")
  
  return(mod_summary)
}

# Model 1: Benthic invertebrate biomass per $cm^2$  
invert_biomass_models <- function(data){
  
  invert_biomass_null_model <- lm(invertebrate_biomass ~ 
                                    1,
                                  data)
  
  null_mod_summary <- extract_model_summary(invert_biomass_null_model, 
                                            "null_mod", 
                                            "invertebrate_biomass")
  
  # version 1
  invert_biomass_model1 <- lm(invertebrate_biomass ~ 
                                total_disturbance +
                                wetted_width +
                                wolmanD50,
                              data)
  
  mod1_summary <- extract_model_summary(invert_biomass_model1, 
                                        "mod1", 
                                        "invertebrate_biomass")
  
  # version 2
  invert_biomass_model2 <- lm(invertebrate_biomass ~ 
                                total_road_density +
                                wetted_width +
                                wolmanD50, 
                              data)  
  mod2_summary <- extract_model_summary(invert_biomass_model2, 
                                        "mod2", 
                                        "invertebrate_biomass")
  
  # version 3
  invert_biomass_model3 <- lm(invertebrate_biomass ~ 
                                high_human_impact + 
                                wetted_width + 
                                wolmanD50, 
                              data)  
  mod3_summary <- extract_model_summary(invert_biomass_model3, 
                                        "mod3", 
                                        "invertebrate_biomass")
  
  # version 4
  invert_biomass_model4 <- lm(invertebrate_biomass ~ 
                                total_road_density + 
                                high_human_impact + 
                                wetted_width + 
                                wolmanD50, 
                              data)  
  mod4_summary <- extract_model_summary(invert_biomass_model4, 
                                        "mod4", 
                                        "invertebrate_biomass")
  
  # version 5
  invert_biomass_model5 <- lm(invertebrate_biomass ~ 
                                total_disturbance + 
                                total_road_density + 
                                wetted_width + 
                                wolmanD50, 
                              data)  
  mod5_summary <- extract_model_summary(invert_biomass_model5, 
                                        "mod5", 
                                        "invertebrate_biomass")
  
  # version 6
  invert_biomass_model6 <- lm(invertebrate_biomass ~ 
                                total_disturbance + 
                                high_human_impact + 
                                wetted_width + 
                                wolmanD50, 
                              data)  
  mod6_summary <- extract_model_summary(invert_biomass_model6, 
                                        "mod6", 
                                        "invertebrate_biomass")
  
  # version 7
  invert_biomass_model7 <- lm(invertebrate_biomass ~ 
                                total_disturbance + 
                                total_road_density + 
                                high_human_impact + 
                                wetted_width + 
                                wolmanD50, 
                              data)  
  mod7_summary <- extract_model_summary(invert_biomass_model7, 
                                        "mod7", 
                                        "invertebrate_biomass")
  
  # version 8 
  invert_biomass_model8 <- lm(invertebrate_biomass ~ 
                                wetted_width + 
                                wolmanD50 + 
                                canopy +
                                flow, 
                              data)
  mod8_summary <- extract_model_summary(invert_biomass_model8, 
                                        "mod8", 
                                        "invertebrate_biomass")
  
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
} 

# Model 2: EPT index  
ept_index_models <- function(data){
  
  # null version
  ept_index_null_model <- lm(ept_index ~ 
                               1, 
                             data)
  
  null_mod_summary <- extract_model_summary(ept_index_null_model, 
                                            "null_mod", 
                                            "ept_index")
  
  # version 1
  ept_index_model1 <- lm(ept_index ~ 
                           total_disturbance + 
                           wetted_width + 
                           wolmanD50, 
                         data)
  
  mod1_summary <- extract_model_summary(ept_index_model1, 
                                        "mod1", 
                                        "ept_index")
  
  # version 2
  ept_index_model2 <- lm(ept_index ~ 
                           total_road_density +
                           wetted_width + 
                           wolmanD50, 
                         data)
  mod2_summary <- extract_model_summary(ept_index_model2, 
                                        "mod2", 
                                        "ept_index")
  
  # version 3
  ept_index_model3 <- lm(ept_index ~ 
                           high_human_impact + 
                           wetted_width + 
                           wolmanD50, 
                         data)
  mod3_summary <- extract_model_summary(ept_index_model3, 
                                        "mod3", 
                                        "ept_index")
  
  # version 4
  ept_index_model4 <- lm(ept_index ~ 
                           total_road_density + 
                           high_human_impact +
                           wetted_width + 
                           wolmanD50, 
                         data)
  mod4_summary <- extract_model_summary(ept_index_model4, 
                                        "mod4", 
                                        "ept_index")
  
  # version 5
  ept_index_model5 <- lm(ept_index ~ 
                           total_disturbance +
                           total_road_density + 
                           wetted_width + 
                           wolmanD50, 
                         data)
  mod5_summary <- extract_model_summary(ept_index_model5, 
                                        "mod5", 
                                        "ept_index")
  
  # version 6
  ept_index_model6 <- lm(ept_index ~ 
                           total_disturbance + 
                           high_human_impact + 
                           wetted_width + 
                           wolmanD50, 
                         data)  
  mod6_summary <- extract_model_summary(ept_index_model6, 
                                        "mod6", 
                                        "ept_index")
  
  # version 7 
  ept_index_model7 <- lm(ept_index ~ 
                           total_disturbance + 
                           total_road_density + 
                           high_human_impact + 
                           wetted_width + 
                           wolmanD50, 
                         data)
  mod7_summary <- extract_model_summary(ept_index_model7, 
                                        "mod7", 
                                        "ept_index")
  
  # version 8 
  ept_index_model8 <- lm(ept_index ~ 
                           wetted_width + 
                           wolmanD50 + 
                           canopy +
                           flow, 
                         data)
  mod8_summary <- extract_model_summary(ept_index_model8, 
                                        "mod8", 
                                        "ept_index")
  
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# Model 3: periphyton biomass
periphyton_models <- function(data){
  
  
  # version null
  periphyton_null_model <- lm(periphyton_biomass ~ 
                                 1, 
                               data)  
  null_mod_summary <- extract_model_summary(periphyton_null_model, 
                                            "null_mod", 
                                            "periphyton_biomass")
  
  # version 1
  periphyton_model1 <- lm(periphyton_biomass ~ 
                             total_disturbance + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod1_summary <- extract_model_summary(periphyton_model1, 
                                        "mod1", 
                                        "periphyton_biomass")
  
  # version 2 (swap depth for wetted width)
  periphyton_model2 <- lm(periphyton_biomass ~ 
                             total_road_density + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod2_summary <- extract_model_summary(periphyton_model2, 
                                        "mod2", 
                                        "periphyton_biomass")
  
  # version 3 (back to wetted width, percent lake)
  periphyton_model3 <- lm(periphyton_biomass ~ 
                             high_human_impact + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod3_summary <- extract_model_summary(periphyton_model3, 
                                        "mod3", 
                                        "periphyton_biomass")
  
  # version 4 (back to percent wetland, add canopy cover)
  periphyton_model4 <- lm(periphyton_biomass ~ 
                             total_road_density + 
                             high_human_impact + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod4_summary <- extract_model_summary(periphyton_model4, 
                                        "mod4", 
                                        "periphyton_biomass")
  
  # version 5
  periphyton_model5 <- lm(periphyton_biomass ~ 
                             total_disturbance +
                             total_road_density + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod5_summary <- extract_model_summary(periphyton_model5, 
                                        "mod5", 
                                        "periphyton_biomass")
  
  # version 6 (back to wetted width, percent lake)
  periphyton_model6 <- lm(periphyton_biomass ~ 
                             total_disturbance + 
                             high_human_impact + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod6_summary <- extract_model_summary(periphyton_model6, 
                                        "mod6", 
                                        "periphyton_biomass")
  
  # version 7 (back to percent wetland, add canopy cover)
  periphyton_model7 <- lm(periphyton_biomass ~ 
                             total_disturbance + 
                             total_road_density + 
                             high_human_impact + 
                             wetted_width + 
                             canopy, 
                           data)  
  mod7_summary <- extract_model_summary(periphyton_model7, 
                                        "mod7", 
                                        "periphyton_biomass")
  
  # version 8
  periphyton_model8 <- lm(periphyton_biomass ~ 
                             wetted_width + 
                             percent_wetland +
                             wolmanD50 + 
                             canopy, 
                           data)  
  mod8_summary <- extract_model_summary(periphyton_model8, 
                                        "mod8", 
                                        "periphyton_biomass")
  
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# Model 4: % shredder 
shredder_models <- function(data){
  
  # version 1
  shredder_null_model <- lm(shredder ~ 
                              1, 
                            data)  
  null_mod_summary <- extract_model_summary(shredder_null_model, 
                                            "null_mod", 
                                            "shredders")
  
  # version 1
  shredder_model1 <- lm(shredder ~ 
                          total_disturbance + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod1_summary <- extract_model_summary(shredder_model1, 
                                        "mod1", 
                                        "shredders")
  
  # version 2
  shredder_model2 <- lm(shredder ~ 
                          total_road_density + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod2_summary <- extract_model_summary(shredder_model2, 
                                        "mod2", 
                                        "shredders")
  
  # version 3
  shredder_model3 <- lm(shredder ~ 
                          high_human_impact + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod3_summary <- extract_model_summary(shredder_model3, 
                                        "mod3", 
                                        "shredders")
  
  # version 4
  shredder_model4 <- lm(shredder ~ 
                          total_road_density + 
                          high_human_impact + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod4_summary <- extract_model_summary(shredder_model4, 
                                        "mod4", 
                                        "shredders")
  
  # version 5
  shredder_model5 <- lm(shredder ~ 
                          total_disturbance + 
                          total_road_density + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod5_summary <- extract_model_summary(shredder_model5, 
                                        "mod5", 
                                        "shredders")
  
  # version 6
  shredder_model6 <- lm(shredder ~ 
                          total_disturbance + 
                          high_human_impact + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod6_summary <- extract_model_summary(shredder_model6, 
                                        "mod6", 
                                        "shredders")
  
  # version 7
  shredder_model7 <- lm(shredder ~ 
                          total_disturbance + 
                          total_road_density + 
                          high_human_impact + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod7_summary <- extract_model_summary(shredder_model7, 
                                        "mod7", 
                                        "shredders")
  
  # version 8
  shredder_model8 <- lm(shredder ~ 
                          canopy + 
                          flow + 
                          wetted_width + 
                          wolmanD50, 
                        data)  
  mod8_summary <- extract_model_summary(shredder_model8, 
                                        "mod8", 
                                        "shredders")
  
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# Model 5: Total nitrogen
nitrogen_models <- function(data){
  
  # version 1
  total_nitrogen_null_model <- lm(tn ~ 
                                    1, 
                                  data)  
  null_mod_summary <- extract_model_summary(total_nitrogen_null_model, 
                                            "null_mod", 
                                            "nitrogen")
  # version 1
  total_nitrogen_model1 <- lm(tn ~ 
                                total_disturbance + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod1_summary <- extract_model_summary(total_nitrogen_model1, 
                                        "mod1", 
                                        "nitrogen")
  
  # version 2
  total_nitrogen_model2 <- lm(tn ~ 
                                total_road_density + 
                                wetted_width + 
                                wetted_width + 
                                percent_wetland,
                              data)  
  mod2_summary <- extract_model_summary(total_nitrogen_model2, 
                                        "mod2", 
                                        "nitrogen")
  
  # version 3
  total_nitrogen_model3 <- lm(tn ~ 
                                high_human_impact + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod3_summary <- extract_model_summary(total_nitrogen_model3, 
                                        "mod3", 
                                        "nitrogen")
  
  # version 4
  total_nitrogen_model4 <- lm(tn ~ 
                                total_road_density + 
                                high_human_impact + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod4_summary <- extract_model_summary(total_nitrogen_model4, 
                                        "mod4", 
                                        "nitrogen")
  # version 5
  total_nitrogen_model5 <- lm(tn ~ 
                                total_disturbance +
                                total_road_density + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod5_summary <- extract_model_summary(total_nitrogen_model5, 
                                        "mod5", 
                                        "nitrogen")
  
  # version 6
  total_nitrogen_model6 <- lm(tn ~ 
                                total_disturbance + 
                                high_human_impact + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod6_summary <- extract_model_summary(total_nitrogen_model6, 
                                        "mod6", 
                                        "nitrogen")
  
  # version 7
  total_nitrogen_model7 <- lm(tn ~ 
                                total_disturbance + 
                                total_road_density + 
                                high_human_impact + 
                                wetted_width + 
                                percent_wetland, 
                              data)  
  mod7_summary <- extract_model_summary(total_nitrogen_model7, 
                                        "mod7", 
                                        "nitrogen")
  # version 8
  total_nitrogen_model8 <- lm(tn ~ 
                                wetted_width + 
                                canopy + 
                                percent_lake +
                                percent_wetland, 
                              data)  
  mod8_summary <- extract_model_summary(total_nitrogen_model8, 
                                        "mod8", 
                                        "nitrogen")
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# Model 6: Specific conductivity
conductivity_models <- function(data){
  
  # version null
  specific_conductivity_null_model <- lm(conductivity ~ 
                                             1, 
                                           data)
  null_mod_summary <- extract_model_summary(specific_conductivity_null_model, 
                                            "null_mod", 
                                            "specific_conductivity")
  
  # version 1
  specific_conductivity_model1 <- lm(conductivity ~ 
                                         total_disturbance + 
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod1_summary <- extract_model_summary(specific_conductivity_model1, 
                                        "mod1", 
                                        "specific_conductivity")
  
  # version 2
  specific_conductivity_model2 <- lm(conductivity ~ 
                                         total_road_density + 
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod2_summary <- extract_model_summary(specific_conductivity_model2, 
                                        "mod2", 
                                        "specific_conductivity")
  
  # version 3
  specific_conductivity_model3 <- lm(conductivity ~ 
                                         high_human_impact +
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod3_summary <- extract_model_summary(specific_conductivity_model3, 
                                        "mod3", 
                                        "specific_conductivity")
  
  # version 4
  specific_conductivity_model4 <- lm(conductivity ~ 
                                         total_road_density +
                                         high_human_impact +
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod4_summary <- extract_model_summary(specific_conductivity_model4, 
                                        "mod4", 
                                        "specific_conductivity")
  
  # version 5
  specific_conductivity_model5 <- lm(conductivity ~ 
                                         total_disturbance + 
                                         total_road_density +
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod5_summary <- extract_model_summary(specific_conductivity_model5, 
                                        "mod5", 
                                        "specific_conductivity")
  
  # version 6
  specific_conductivity_model6 <- lm(conductivity ~ 
                                         total_disturbance + 
                                         high_human_impact + 
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod6_summary <- extract_model_summary(specific_conductivity_model6, 
                                        "mod6", 
                                        "specific_conductivity")
  
  # version 7
  specific_conductivity_model7 <- lm(conductivity ~ 
                                         total_disturbance + 
                                         total_road_density + 
                                         high_human_impact + 
                                         wetted_width + 
                                         percent_wetland, 
                                       data)
  mod7_summary <- extract_model_summary(specific_conductivity_model7, 
                                        "mod7", 
                                        "specific_conductivity")
  
  # version 8
  specific_conductivity_model8 <- lm(conductivity ~ 
                                         wetted_width + 
                                         canopy + 
                                         percent_wetland + 
                                         percent_lake, 
                                       data)
  mod8_summary <- extract_model_summary(specific_conductivity_model8, 
                                        "mod8", 
                                        "specific_conductivity")
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# Model 7: Substrate embeddedness
embeddedness_models <- function(data){
  
  # version 1
  percent_embeddedness_null_model <- lm(embeddedness ~ 
                                         1, 
                                       data)  
  null_mod_summary <- extract_model_summary(percent_embeddedness_null_model, 
                                            "null_mod", 
                                            "embeddedness")
  
  # version 1
  percent_embeddedness_model1 <- lm(embeddedness ~ 
                                     total_disturbance + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod1_summary <- extract_model_summary(percent_embeddedness_model1, 
                                        "mod1", 
                                        "embeddedness")
  
  # version 2
  percent_embeddedness_model2 <- lm(embeddedness ~ 
                                     total_road_density + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod2_summary <- extract_model_summary(percent_embeddedness_model2, 
                                        "mod2", 
                                        "embeddedness")
  
  # version 3
  percent_embeddedness_model3 <- lm(embeddedness ~ 
                                     high_human_impact + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod3_summary <- extract_model_summary(percent_embeddedness_model3, 
                                        "mod3", 
                                        "embeddedness")
  
  # version 4
  percent_embeddedness_model4 <- lm(embeddedness ~ 
                                     total_road_density + 
                                     high_human_impact + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod4_summary <- extract_model_summary(percent_embeddedness_model4, 
                                        "mod4", 
                                        "embeddedness")
  
  # version 5
  percent_embeddedness_model5 <- lm(embeddedness ~ 
                                     total_disturbance + 
                                     total_road_density + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod5_summary <- extract_model_summary(percent_embeddedness_model5, 
                                        "mod5", 
                                        "embeddedness")
  
  # version 6
  percent_embeddedness_model6 <- lm(embeddedness ~ 
                                     total_disturbance + 
                                     high_human_impact + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod6_summary <- extract_model_summary(percent_embeddedness_model6, 
                                        "mod6", 
                                        "embeddedness")
  
  # version 7
  percent_embeddedness_model7 <- lm(embeddedness ~ 
                                     total_disturbance + 
                                     total_road_density + 
                                     high_human_impact + 
                                     wolmanD50 + 
                                     percent_lake, 
                                   data)  
  mod7_summary <- extract_model_summary(percent_embeddedness_model7, 
                                        "mod7", 
                                        "embeddedness")
  
  # version 8
  percent_embeddedness_model8 <- lm(embeddedness ~ 
                                     wetted_width + 
                                     wolmanD50 + 
                                     percent_lake +
                                     flow, 
                                   data)  
  mod8_summary <- extract_model_summary(percent_embeddedness_model8, 
                                        "mod8", 
                                        "embeddedness")
  
  # combine dataframes for all models
  model_summary_df <- bind_rows(null_mod_summary,
                                mod1_summary,
                                mod2_summary,
                                mod3_summary,
                                mod4_summary,
                                mod5_summary,
                                mod6_summary,
                                mod7_summary,
                                mod8_summary)
  
  return(model_summary_df)
}

# 4) Compare all models at each spatial extent ----

# "all models" function
run_all_models <- function(data){
  
  # run all models and convert output to data frame
  invert_model_output <- invert_biomass_models(data)
  invert_model_df <- invert_model_output[1:10] |> data.frame()
  
  ept_model_output <- ept_index_models(data)
  ept_model_df <- ept_model_output[1:10] |> data.frame()
  
  periphyton_model_output <- periphyton_models(data)
  periphyton_model_df <- periphyton_model_output[1:10] |> data.frame()
  
  shredder_model_output <- shredder_models(data)
  shredder_model_df <- shredder_model_output[1:10] |> data.frame()
  
  nitrogen_model_output <- nitrogen_models(data)
  nitrogen_model_df <- nitrogen_model_output[1:10] |> data.frame() 
  
  conductivity_model_output <- conductivity_models(data)
  conductivity_model_df <- conductivity_model_output[1:10] |> data.frame() 
  
  embeddedness_model_output <- embeddedness_models(data)
  embeddedness_model_df <- embeddedness_model_output[1:10] |> data.frame() 
  
  # combine all data frames together
  all_models_df <- bind_rows(invert_model_df, 
                             ept_model_df, 
                             periphyton_model_df, 
                             shredder_model_df,
                             nitrogen_model_df, 
                             conductivity_model_df,
                             embeddedness_model_df)
  
  return(all_models_df)
}

model_summary_large <- run_all_models(rescaled_data_large)
model_summary_med <- run_all_models(rescaled_data_med)
model_summary_small <- run_all_models(rescaled_data_small)

# first combine data for all scales
model_summary_large$scale = "catchment"
model_summary_med$scale = "riparian"
model_summary_small$scale = "local"

model_summary <- rbind(model_summary_large, 
                       model_summary_med, 
                       model_summary_small)

# export model summary

model_data_filtered <-model_summary %>% 
  filter(scale=="catchment") %>%
  group_by(model) %>% 
  mutate(delta_aic = aic - min(aic)) %>%
  filter(delta_aic <= 2)

write.csv(model_data_filtered, "output/best_catchmentv3_models_.csv", row.names=FALSE)
write.csv(model_summary, "output/all_model_resultsv3.csv", row.names=FALSE)

# 5) Format and export model summary ----
# stats  model summary
library(dplyr)
library(stringr)
library(reshape2)

# read in csv with model results
setwd("C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc")
model_summary <- read.csv("output/all_model_resultsv3.csv")

# make new column combining estimate and standard deviation with a ± between
formatted_model_output <- model_summary %>% 
  mutate(full_estimate = paste(as.character(round(estimate, 2)), "(", 
                               as.character(round(std_error, 2)), ")")) %>%
  mutate(across("predictor_var", str_replace, "total_disturbance", "forest disturbance")) %>%
  mutate(across("predictor_var", str_replace, "total_road_density", "road density")) %>% 
  mutate(across("predictor_var", str_replace, "high_human_impact", "human impact")) %>% 
  mutate(across("predictor_var", str_replace, "wetted_width", "wetted width")) %>% 
  mutate(across("predictor_var", str_replace, "wolmanD50", "substrate size")) %>%
  mutate(across("predictor_var", str_replace, "percent_lake", "% lake")) %>% 
  mutate(across("predictor_var", str_replace, "percent_wetland", "% wetland")) %>% 
  mutate(across("predictor_var", str_replace, "percent_barren", "% barren")) %>%
  mutate(across("predictor_var", str_replace, "(Intercept)", "intercept")) %>% 
  mutate(across("model", str_replace, "periphyton_biomass", "periphyton biomass")) %>% 
  mutate(across("model", str_replace, "ept_index", "EPT index")) %>% 
  mutate(across("model", str_replace, "invertebrate_biomass", "invertebrate biomass")) %>% 
  mutate(across("model", str_replace, "nitrogen", "total nitrogen")) %>% 
  mutate(across("model", str_replace, "shredders", "% shredders")) %>% 
  mutate(across("scale", str_replace, "sub-catchment", "catchment")) %>%
  mutate(across("model_num", str_replace, "null_mod", "null")) %>%
  mutate(across("model_num", str_replace, "mod1", "1")) %>%
  mutate(across("model_num", str_replace, "mod2", "2")) %>%
  mutate(across("model_num", str_replace, "mod3", "3")) %>%
  mutate(across("model_num", str_replace, "mod4", "4")) %>%
  mutate(across("model_num", str_replace, "mod5", "5")) %>%
  mutate(across("model_num", str_replace, "mod6", "6")) %>%
  mutate(across("model_num", str_replace, "mod7", "7")) %>%
  mutate(across("model_num", str_replace, "mod8", "8")) %>%
  group_by(scale, model) %>% 
  mutate(delta_aic = aic - min(aic)) %>% 
  select(scale, 
         model, 
         model_num, 
         predictor_var,
         full_estimate, 
         log_like, 
         delta_aic,
         aic,
         adj_r2) %>% 
  rename("estimate" = "full_estimate", 
         "model number" = "model_num", 
         "predictor" = "predictor_var", 
         "log likelihood" = "log_like",
         "extent" = "scale")

write.csv(formatted_model_output, "output/formatted_stats_model_output_v3.csv", row.names=FALSE)

# 6) Plot results from the top models ----

# first replace all column names with a clean label
plotting_data <- model_summary %>% 
  mutate(dplyr::across("predictor_var", str_replace, "total_disturbance", "forest disturbance")) %>% 
  mutate(dplyr::across("predictor_var", str_replace, "total_road_density", "unpaved road density")) %>% 
  mutate(dplyr::across("predictor_var", str_replace, "high_human_impact", "human impact index")) %>% 
  mutate(dplyr::across("predictor_var", str_replace, "wetted_width", "wetted width")) %>% 
  mutate(dplyr::across("predictor_var", str_replace, "wolmanD50", "substrate size")) %>%
  mutate(dplyr::across("predictor_var", str_replace, "percent_lake", "% lake")) %>% 
  mutate(dplyr::across("predictor_var", str_replace, "percent_wetland", "% wetland")) %>% 
  mutate(dplyr::across("model", str_replace, "periphyton_biomass", "periphyton biomass")) %>% 
  mutate(dplyr::across("model", str_replace, "specific_conductivity", "specific conductivity")) %>% 
  mutate(dplyr::across("model", str_replace, "ept_index", "EPT index")) %>% 
  mutate(dplyr::across("model", str_replace, "invertebrate_biomass", "invertebrate biomass")) %>% 
  mutate(dplyr::across("model", str_replace, "nitrogen", "total dissolved nitrogen")) %>% 
  mutate(dplyr::across("model", str_replace, "shredders", "% shredders"))

# plot on vertical histogram
# 1) all spatial extents together
invert_bar_data <- plotting_data %>% 
  # filter(model_num != "null_mod") %>%
  group_by(model, scale) %>% # 
  mutate(delta_aic = aic - min(aic)) %>%
  filter(aic == min(aic)) %>%
  # filter(delta_aic <= 2) %>%
  filter(predictor_var != "(Intercept)")
invert_bar_data$significant <- dplyr::case_when(
  invert_bar_data$p_val < 0.05 ~ TRUE,
  invert_bar_data$p_val >= 0.05 ~ FALSE
)

bar_data_local <- invert_bar_data %>% filter(scale=="local")
bar_data_riparian <- invert_bar_data %>% filter(scale=="riparian") %>% filter(model != "invertebrate biomass")
bar_data_catchment <- invert_bar_data %>% filter(scale=="catchment") %>% filter(model != "invertebrate biomass")

model_summary_plot <- ggplot(bar_data_local, 
                             aes(predictor_var, 
                                 estimate, 
                                 fill=predictor_var)) +
  geom_col()+
  # scale_y_reverse(limits=c(1,-1)) + 
  theme_classic() + 
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) + 
  coord_flip() + 
  geom_hline(yintercept=0, 
             linetype="solid", 
             color="black", 
             size=0.5) + 
  ylab("estimated slope") + 
  xlab("predictor") + 
  scale_fill_viridis_d() + 
  facet_grid(cols = vars(model)) +  # , rows =vars(scale)
  theme(strip.text = element_text(
    size = 8)) + 
  geom_text(aes(label = ifelse(significant, "*", "")),  # , group = model
            size = 10 / .pt) + 
  geom_errorbar(aes(ymin=estimate-std_error, ymax=estimate+std_error), width=.2,
                position=position_dodge(.9))

# export
ggsave(file="C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc/output/model_summary_local_v3.svg", plot=model_summary_plot, width=5, height=2)
