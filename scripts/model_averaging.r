# This script is used to perform model averaging for all models in the
# stats models script after screening for uninformative parameters.

# Sections:
# 1) Read in rescaled data for each spatial extent
# 2) Create functions to extract data from model summaries
# 3) Create functions for each model comparison, including only informative models
# 4) Compare models at each spatial extent
# 5) Format and export results
# 6) Generate plots

# 1) Read in rescaled data ----

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tibble)
library(stringr)
library(wiqid)
library(nlme)
library(AICcmodavg)

# read in rescaled data
rescaled_data_local <- read.csv("output/rescaled_data_local.csv")
rescaled_data_riparian <- read.csv("output/rescaled_data_riparian.csv")
rescaled_data_catchment <- read.csv("output/rescaled_data_catchment.csv")

# 2) Create functions to extract data from model summaries ----

# function to extract data from the model average summary
extract_model_average_data <- function(model_average_output, model_name, extent, param_name) {

  #'@description this function takes the output from the modavg function
  #'and returns it as a formatted dataframe
  #'@param model_average_output output from modavg function
  #'@param model_name name of the response vartiable in the model set,
  #' can be "ept_index", "invertebrate_biomass", "periphyton_biomass",
  #' "specific_conductivity", "nitrogen", "embeddedness", or "shredders"
  #'@param extent spatial extent of the data used to create
  #' model_average_output, can be "local", "riparian", or "catchment"
  #'@param param_name parameter used to generate model_average_output
  #'@return a formatted data frame with key components of the model summary

  # extract items from the model output
  extent <- extent
  model <- model_name
  parameter <- param_name
  beta_est <- model_average_output$Mod.avg.beta
  st_error <- model_average_output$Uncond.SE
  conf_level <- model_average_output$Conf.level
  lower_cl <- model_average_output$Lower.CL
  upper_cl <- model_average_output$Upper.CL

  data_table <- data.frame(extent,
                           model,
                           parameter,
                           beta_est,
                           st_error,
                           conf_level,
                           lower_cl,
                           upper_cl)

  return(data_table)
}

# function to extract data from the null models
extract_null_summary <- function(null_model_output, extent, model_name) {

  #'@description this function takes the output from the lm function
  #'and returns it as a formatted dataframe
  #'@param null_model_output output from lm function
  #'@param extent spatial extent of the data used to create the
  #' null model output, can be "local", "riparian", or "catchment"
  #'@param model_name name of the response vartiable in the model set,
  #' can be "ept_index", "invertebrate_biomass", "periphyton_biomass",
  #' "specific_conductivity", "nitrogen", "embeddedness", or "shredders"
  #'@return a formatted data frame with key components of the model output

  # extract items from the model output
  model_summary <- summary(null_model_output)

  # select coefficients and convert to dataframe
  model_summary_df <- model_summary$coefficients %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "predictor_var") %>%
    mutate(extent = extent,
           model = model_name,
           parameter = "intercept",
           conf_level = NA,
           lower_cl = NA,
           upper_cl = NA) %>%
    rename("beta_est" = "Estimate",
           "st_error" = "Std..Error") %>%
    select(extent,
           model,
           parameter,
           beta_est,
           st_error,
           conf_level,
           lower_cl,
           upper_cl)

  return(model_summary_df)
}

# 3) Create functions for each model comparison ----

## Model 1: Benthic invertebrate biomass per cm^2 ----
invert_biomass_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with benthic invertebrate biomass as the response variable, and
  #' returns a formatted dataframe of the results by calling on the "extract
  #' model average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local") {
    invert_biomass_null_model <- lm(invertebrate_biomass ~
                                      1,
                                    data)

    invert_biomass_model10 <- lm(invertebrate_biomass ~
                                   total_road_density,
                                 data)

    # make a list of model outputs
    model_list <- list("null_mod" = invert_biomass_null_model,
                       "mod10" = invert_biomass_model10)

    # compute model averaging
    param_name <-  "total_road_density"
    averaged_model <- modavg(cand.set = model_list, parm = param_name)
    model_table <- extract_model_average_data(averaged_model,
                                              "invertebrate_biomass",
                                              extent,
                                              param_name)
  }

  if (extent == "riparian" || extent == "catchment") {
    invert_biomass_null_model <- lm(invertebrate_biomass ~
                                      1,
                                    data)

    # make a list of model outputs
    model_list <- invert_biomass_null_model

    # extract data from the null model
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "invertebrate_biomass")
  }

  return(model_table)
}

## Model 2: EPT index ----
ept_index_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with EPT index as the response variable, and returns a formatted
  #' dataframe of the results by calling on the "extract model average data"
  #' and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local") {

    ept_index_null_model <- lm(ept_index ~
                                 1,
                               data)

    ept_index_model9 <- lm(ept_index ~
                             total_disturbance,
                           data)

    model_list <- list("null_mod" = ept_index_null_model,
                       "mod9" = ept_index_model9)

    # compute model averaging
    param_name <- "total_disturbance"
    averaged_model <- modavg(cand.set = model_list, parm = param_name)
    model_table <- extract_model_average_data(averaged_model,
                                              "ept_index",
                                              extent,
                                              param_name)
  }

  if (extent == "riparian") {

    ept_index_null_model <- lm(ept_index ~
                                 1,
                               data)

    ept_index_model9 <- lm(ept_index ~
                             total_disturbance,
                           data)

    ept_index_model10 <- lm(ept_index ~
                              total_road_density,
                            data)

    model_list <- list("null_mod" = ept_index_null_model,
                       "mod9" = ept_index_model9,
                       "mod10" = ept_index_model10)

    # compute model averaging
    param_name1 <- "total_road_density"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "ept_index",
                                               extent,
                                               param_name1)

    param_name2 <- "total_disturbance"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "ept_index",
                                               extent,
                                               param_name2)

    model_table <- bind_rows(model_table1,
                             model_table2)
  }

  if (extent == "catchment") {

    ept_index_null_model <- lm(ept_index ~
                                 1,
                               data)

    model_list <- ept_index_null_model

    # extract data from the null model
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "ept_index")
  }

  return(model_table)
}

## Model 3: periphyton ----
periphyton_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with periphyton biomass as the response variable, and returns a
  #' formatted dataframe of the results by calling on the "extract model
  #' average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local" || extent == "riparian" || extent == "catchment") {

    periphyton_null_model <- lm(periphyton_biomass ~
                                  1,
                                data)

    model_list <- periphyton_null_model

    # extract data from the null model
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "periphyton_biomass")
  }

  return(model_table)
}

## Model 4: % shredder ----
shredder_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with percent shredders as the response variable, and returns a
  #' formatted dataframe of the results by calling on the "extract model
  #' average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local") {

    shredder_null_model <- lm(shredder ~
                                1,
                              data)

    shredder_model1 <- lm(shredder ~
                            total_disturbance +
                              wetted_width +
                              wolmanD50,
                          data)

    shredder_model3 <- lm(shredder ~
                            high_human_impact +
                              wetted_width +
                              wolmanD50,
                          data)

    shredder_model4 <- lm(shredder ~
                            total_road_density +
                              high_human_impact +
                              wetted_width +
                              wolmanD50,
                          data)

    shredder_model8 <- lm(shredder ~
                            wetted_width +
                              wolmanD50,
                          data)

    shredder_model9 <- lm(shredder ~
                            total_disturbance,
                          data)

    shredder_model11 <- lm(shredder ~
                             high_human_impact,
                           data)

    model_list <- list("null_mod" = shredder_null_model,
                       "mod1" = shredder_model1,
                       "mod3" = shredder_model3,
                       "mod4" = shredder_model4,
                       "mod8" = shredder_model8,
                       "mod9" = shredder_model9,
                       "mod11" = shredder_model11)

    # compute model averaging
    param_name1 <- "total_road_density"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "shredders",
                                               extent,
                                               param_name1)

    param_name2 <- "total_disturbance"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "shredders",
                                               extent,
                                               param_name2)

    param_name3 <- "wetted_width"
    averaged_model3 <- modavg(cand.set = model_list, parm = param_name3)
    model_table3 <- extract_model_average_data(averaged_model3,
                                               "shredders",
                                               extent,
                                               param_name3)

    param_name4 <- "wolmanD50"
    averaged_model4 <- modavg(cand.set = model_list, parm = param_name4)
    model_table4 <- extract_model_average_data(averaged_model4,
                                               "shredders",
                                               extent,
                                               param_name4)

    param_name5 <- "high_human_impact"
    averaged_model5 <- modavg(cand.set = model_list, parm = param_name5)
    model_table5 <- extract_model_average_data(averaged_model5,
                                               "shredders",
                                               extent,
                                               param_name5)

    model_table <- bind_rows(model_table1,
                             model_table2,
                             model_table3,
                             model_table4,
                             model_table5)
  }

  if (extent == "riparian") {

    shredder_null_model <- lm(shredder ~ 1, data)

    shredder_model1 <- lm(shredder ~
                            total_disturbance +
                              wetted_width +
                              wolmanD50,
                          data)

    shredder_model2 <- lm(shredder ~
                            total_road_density +
                              wetted_width +
                              wolmanD50,
                          data)

    shredder_model8 <- lm(shredder ~
                            wetted_width +
                              wolmanD50,
                          data)

    shredder_model9 <- lm(shredder ~
                            total_disturbance,
                          data)

    model_list <- list("null_mod" = shredder_null_model,
                       "mod1" = shredder_model1,
                       "mod2" = shredder_model2,
                       "mod8" = shredder_model8,
                       "mod9" = shredder_model9)

    # compute model averaging
    param_name1 <- "total_road_density"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "shredders",
                                               extent,
                                               param_name1)

    param_name2 <- "total_disturbance"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "shredders",
                                               extent,
                                               param_name2)

    param_name3 <- "wetted_width"
    averaged_model3 <- modavg(cand.set = model_list, parm = param_name3)
    model_table3 <- extract_model_average_data(averaged_model3,
                                               "shredders",
                                               extent,
                                               param_name3)

    param_name4 <- "wolmanD50"
    averaged_model4 <- modavg(cand.set = model_list, parm = param_name4)
    model_table4 <- extract_model_average_data(averaged_model4,
                                               "shredders",
                                               extent,
                                               param_name4)

    model_table <- bind_rows(model_table1,
                             model_table2,
                             model_table3,
                             model_table4)
  }

  if (extent == "catchment") {

    shredder_null_model <- lm(shredder ~
                                1,
                              data)

    shredder_model8 <- lm(shredder ~
                            wetted_width +
                              wolmanD50,
                          data)

    model_list <- list("null_mod" = shredder_null_model,
                       "mod8" = shredder_model8)

    # compute model averaging
    param_name1 <- "wetted_width"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "shredders",
                                               extent,
                                               param_name1)

    param_name2 <- "wolmanD50"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "shredders",
                                               extent,
                                               param_name2)

    model_table <- bind_rows(model_table1,
                             model_table2)
  }

  return(model_table)
}

## Model 5: Total nitrogen ----
nitrogen_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with total nitrogen as the response variable, and returns a
  #' formatted dataframe of the results by calling on the "extract model
  #' average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local" || extent == "riparian" || extent == "catchment") {

    total_nitrogen_null_model <- lm(tn ~
                                      1,
                                    data)

    model_list <- total_nitrogen_null_model

    # extract data from the null model
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "nitrogen")
  }

  return(model_table)
}

## 6 Specific conductivity ----
conductivity_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with specific conductivity as the response variable, and returns a
  #' formatted dataframe of the results by calling on the "extract model
  #' average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local") {

    specific_conductivity_null_model <- lm(conductivity ~
                                             1,
                                           data)

    model_list <- specific_conductivity_null_model

    # extract data from the null model
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "specific_conductivity")
  }

  if (extent == "riparian") {

    specific_conductivity_null_model <- lm(conductivity ~
                                             1,
                                           data)

    specific_conductivity_model9 <- lm(conductivity ~
                                         total_disturbance,
                                       data)

    specific_conductivity_model11 <- lm(conductivity ~
                                          high_human_impact,
                                        data)

    model_list <- list("null_mod" = specific_conductivity_null_model,
                       "mod9" = specific_conductivity_model9,
                       "mod11" = specific_conductivity_model11)

    # compute model averaging
    param_name1 <- "total_disturbance"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "specific_conductivity",
                                               extent,
                                               param_name1)

    param_name2 <- "high_human_impact"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "specific_conductivity",
                                               extent,
                                               param_name2)

    model_table <- bind_rows(model_table1,
                             model_table2)
  }

  if (extent == "catchment") {

    specific_conductivity_null_model <- lm(conductivity ~
                                             1,
                                           data)

    specific_conductivity_model10 <- lm(conductivity ~
                                          total_road_density,
                                        data)

    specific_conductivity_model11 <- lm(conductivity ~
                                          high_human_impact,
                                        data)

    model_list <- list("null_mod" = specific_conductivity_null_model,
                       "mod10" = specific_conductivity_model10,
                       "mod11" = specific_conductivity_model11)

    # compute model averaging
    param_name1 <- "total_road_density"
    averaged_model1 <- modavg(cand.set = model_list, parm = param_name1)
    model_table1 <- extract_model_average_data(averaged_model1,
                                               "specific_conductivity",
                                               extent,
                                               param_name1)

    param_name2 <- "high_human_impact"
    averaged_model2 <- modavg(cand.set = model_list, parm = param_name2)
    model_table2 <- extract_model_average_data(averaged_model2,
                                               "specific_conductivity",
                                               extent,
                                               param_name2)

    model_table <- bind_rows(model_table1,
                             model_table2)
  }

  return(model_table)
}

## Model 7: Substrate embeddedness ----
embeddedness_models <- function(data, extent) {

  #'@description this function takes rescaled data from one of three spatial
  #' extents, performs model averaging on a set of general linear models
  #' with substrate embeddedness as the response variable, and returns a
  #' formatted dataframe of the results by calling on the "extract model
  #' average data" and "extract null summary" functions
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame with key components of the model output

  if (extent == "local" || extent == "riparian" || extent == "catchment") {

    percent_embeddedness_null_model <- lm(embeddedness ~
                                            1,
                                          data)

    model_list <- percent_embeddedness_null_model

    # extract null summary
    model_table <- extract_null_summary(model_list,
                                        extent,
                                        "embeddedness")
  }

  return(model_table)
}

# 4) Compare models at each spatial extent ----

# make a function to run all models
run_all_models <- function(data, extent) {

  #'@description this function runs the functions for the seven model sets
  #' in the code above
  #'@param data rescaled data at the specified spatial extent, used as
  #' variables in the general linear models
  #'@param extent spatial extent of the data used in the models,
  #' can be "local", "riparian", or "catchment"
  #'@return a formatted data frame compiled from all seven model sets

  # run all models
  invert_models <- invert_biomass_models(data, extent)
  ept_models <- ept_index_models(data, extent)
  periphyton_models <- periphyton_models(data, extent)
  shredder_models <- shredder_models(data, extent)
  nitrogen_models <- nitrogen_models(data, extent)
  conductivity_models <- conductivity_models(data, extent)
  embeddedness_models <- embeddedness_models(data, extent)

  # compile model outputs
  compiled_model_output <- bind_rows(invert_models,
                                     ept_models,
                                     periphyton_models,
                                     shredder_models,
                                     nitrogen_models,
                                     conductivity_models,
                                     embeddedness_models)

  return(compiled_model_output)
}

# run models for each spatial extent
model_summary_catchment <- run_all_models(rescaled_data_catchment, "catchment")
model_summary_riparian <- run_all_models(rescaled_data_riparian, "riparian")
model_summary_local <- run_all_models(rescaled_data_local, "local")

# combine results into one dataframe
all_extent_summary <- bind_rows(model_summary_local,
                                model_summary_riparian,
                                model_summary_catchment)

# 5) Format and export results ----
formatted_model_output <- all_extent_summary %>%
  mutate(parameter = case_when(parameter %in% "total_disturbance" ~ "forest disturbance",
                               parameter %in% "total_road_density" ~ "road density",
                               parameter %in% "high_human_impact" ~ "human impact index",
                               parameter %in% "wetted_width" ~ "wetted width",
                               parameter %in% "wolmanD50" ~ "substrate size",
                               parameter %in% "percent_lake" ~ "% lake",
                               parameter %in% "percent_wetland" ~ "% wetland",
                               parameter %in% "percent_barren" ~ "% barren",
                               parameter %in% "(Intercept)" ~ "intercept")) %>%
  mutate(model =  case_when(model %in% "periphyton_biomass" ~ "periphyton biomass",
                            model %in% "ept_index" ~ "EPT index",
                            model %in% "invertebrate_biomass" ~ "invertebrate biomass",
                            model %in% "specific_conductivity" ~ "specific conductivity",
                            model %in% "nitrogen" ~ "total nitrogen",
                            model %in% "embeddedness" ~ "embeddedness",
                            model %in% "shredders" ~ "% shredders")) %>%
  select(extent,
         model,
         parameter,
         beta_est,
         st_error,
         conf_level,
         lower_cl,
         upper_cl)

# export the formatted results
write.csv(formatted_model_output,
          "output/model_averaging_formatted.csv",
          row.names = FALSE)

# 6) Generate plots ----
plotting_data <- formatted_model_output %>%
  filter(parameter != "intercept")

# separate into a new dataframe fom each extent
plotting_data_local <- plotting_data %>% filter(extent == "local")
plotting_data_riparian <- plotting_data %>% filter(extent == "riparian")
plotting_data_catchment <- plotting_data %>% filter(extent == "catchment")

# generate a plot for each spatial extent
p1 <- ggplot(plotting_data_local, aes(parameter,
                                      beta_est,
                                      fill = parameter)) +
  geom_col() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  coord_flip() +
  geom_hline(yintercept = 0,
             linetype = "solid",
             color = "black",
             linewidth = 0.5) +
  ylab("") +
  xlab("") +
  scale_fill_viridis_d() +
  facet_grid(cols = vars(model)) +
  theme(strip.text = element_text(size = 8)) +
  geom_errorbar(aes(ymin = beta_est - st_error,
                    ymax = beta_est + st_error),
                width = .2,
                position = position_dodge(.9))

p2 <- ggplot(plotting_data_riparian, aes(parameter,
                                         beta_est,
                                         fill = parameter)) +
  geom_col() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  coord_flip() +
  geom_hline(yintercept = 0,
             linetype = "solid",
             color = "black",
             linewidth = 0.5) +
  ylab("") +
  xlab("") +
  scale_fill_viridis_d() +
  facet_grid(cols = vars(model)) +
  theme(strip.text = element_text(size = 8)) +
  geom_errorbar(aes(ymin = beta_est - st_error,
                    ymax = beta_est + st_error),
                width = .2,
                position = position_dodge(.9))

p3 <- ggplot(plotting_data_catchment, aes(parameter,
                                          beta_est,
                                          fill = parameter)) +
  geom_col() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  coord_flip() +
  geom_hline(yintercept = 0,
             linetype = "solid",
             color = "black",
             linewidth = 0.5) +
  ylab("") +
  xlab("") +
  scale_fill_viridis_d() +
  facet_grid(cols = vars(model)) +
  theme(strip.text = element_text(size = 8)) +
  geom_errorbar(aes(ymin = beta_est - st_error,
                    ymax = beta_est + st_error),
                width = .2,
                position = position_dodge(.9))

# export plots
ggsave(file = "output/model_summary_local_averaged.svg",
       plot = p1,
       width = 8,
       height = 2)
ggsave(file = "output/model_summary_riparian_averaged.svg",
       plot = p2,
       width = 8,
       height = 2)
ggsave(file = "output/model_summary_catchment_averaged.svg",
       plot = p3,
       width = 6,
       height = 1.75)
