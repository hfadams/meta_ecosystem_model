# Generate correlogram of empirical data before developing general linear 
# models for statistical analysis

# Load packages
library(GGally)
library(dplyr)
library(ggplot2)
library(svglite)

# read in data
disturbance_data <- read.csv("output/disturbance_data_catchment.csv") %>% 
  dplyr::select(site,
                total_disturbance,
                unpaved_road_density,
                high_human_impact,
                percent_barren,
                percent_lake,
                percent_wetland) %>%
  dplyr::rename(forest_disturbance = total_disturbance)

field_data <- read.csv("output/empirical_stream_data.csv") %>% 
  dplyr::select(site,
                flow,
                wetted_width,
                wolmanD50,
                depth,
                canopy) %>%
  dplyr::rename(substrate_size = wolmanD50)

# combine into one dataset
site_data <- field_data %>% 
  full_join(disturbance_data, by="site") %>% 
  dplyr::select(-site) 

# generate correlogram
correlogram <- ggpairs(site_data)
ggsave(file="output/disturbance_correlogram.svg", plot=correlogram, width=16, height=16)
