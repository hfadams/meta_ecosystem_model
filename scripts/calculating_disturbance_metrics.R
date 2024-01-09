# script for calculating metrics for disturbance and land cover within three spatial extents around each study site. 
# Land cover data are not provided as per our data sharing agreement with the Government of Newfoundland and Parks Canada.

# sections:
# 1) read in shapefiles
# 2) calculate disturbance and land use in each catchment
# 3) code for generating smallest site buffer

# load packages
library(raster)
library(rgdal)
library(sf)
library(elevatr)
library(dplyr)
library(lwgeom)
library(reshape2)
library(devtools)

# 1) Read in shapefiles ----
# transforming all shapefiles to match the stream lengths shapefile crs

# generated in QGIS
stream_reach <- st_read("stream_reach.shp")
catchments <- st_transform(st_read("catchments.shp"), crs(stream_reach)) %>%
  arrange(site)
stream_sites <- st_transform(st_read("stream_sites.shp"), crs(stream_reach))
# drop geometry to select dataframe with site names
site_list <- st_drop_geometry((stream_sites %>% dplyr::select(park, site)), geometry=NULL) %>% 
  distinct()
# read in 100 m buffer around the streams (medium spatial scale)
riparian_extent <- st_transform(st_read("riparian_extent.shp"), crs(stream_reach)) %>% 
  st_zm(drop = TRUE, what = "ZM")
# read in buffer for closest radius from sampling site that captures 10% of the catchment (smallest spatial scale)
site_extent <- st_transform(st_read("local_extent.shp"), crs(stream_reach)) %>% 
  st_zm(drop = TRUE, what = "ZM") %>% rename(catchment_area = ctchmn_) # see code for generating this at the bottom of the script

# digitized by hand
paved_roads <- st_transform(st_read("paved_roads.shp"), crs(stream_reach))
unpaved_roads <- st_transform(st_read("unpaved_roads.shp"), crs(stream_reach))
trails <- st_transform(st_read("trails.shp"), crs(stream_reach))
forest_disturbance <- st_transform(st_read("forest_disturbance.shp"), crs(stream_reach)) %>% 
  st_zm(drop = TRUE, what = "ZM")
forest_disturbance <- st_make_valid(forest_disturbance)
# separate out each type of disturbance from moose mediated meadow layer
insect <- subset(forest_disturbance, disturbanc == "insect")
logging <- subset(forest_disturbance, disturbanc == "logging")
cleared <- subset(forest_disturbance, disturbanc == "cleared")
fire <- subset(forest_disturbance, disturbanc == "fire")

# FRI data (generated from compiling and editing parks and provincial shapefiles)
lakes <- st_transform(st_read("lakes.shp"), crs(stream_reach)) %>% 
  dplyr::select(geometry) %>% 
  st_zm(drop = TRUE, what = "ZM")
barrens <- st_transform(st_read("barrens.shp"), crs(stream_reach)) %>% 
  dplyr::select(geometry) %>% 
  st_zm(drop = TRUE, what = "ZM")
wetland <- st_transform(st_read("wetland.shp"), crs(stream_reach)) %>% 
  dplyr::select(geometry) %>% 
  st_zm(drop = TRUE, what = "ZM")

# read in human footprint (converted to shapefile using "vectorize" tool)
human_footprint <- st_transform(st_read("human_footrpint.shp"), crs(stream_reach))

# 2) calculate disturbance and land use in each catchment ----
# calculate area for each catchment
catchments$catchment_area <- st_area(catchments$geometry)
riparian_extent$catchment_area <- st_area(riparian_extent$geometry)
site_extent$catchment_area <- st_area(site_extent$geometry)

# flatten polygons to be 2D
sf_use_s2(FALSE)

# function to calculate % area of polygons
percent_of_catchment <- function(shapefile, spatial_extent_shapefile){
  # clip to intersection of spatial extent shapefile
  clipped_layer <- st_intersection(shapefile, spatial_extent_shapefile)
  # calculate area of each polygon
  clipped_layer$shapefile_area <- st_area(clipped_layer$geometry)
  # convert from spatial polygon dataframe to dataframe without geometry
  shapefile_df <- as.data.frame(clipped_layer)
  # group by site and calculate percent area covered
  shapefile_percent <- shapefile_df %>% 
    group_by(site) %>% 
    mutate(percent_area = sum(shapefile_area)/catchment_area*100) %>% 
    dplyr::select(site, catchment_area, percent_area) %>% 
    distinct()
  
  return(shapefile_percent) 
}

# function to calculate total road/trail length
density_in_catchment <- function(shapefile, spatial_extent_shapefile){
  # clip to intersection of spatial extent shapefile
  clipped_layer <- st_intersection(shapefile, spatial_extent_shapefile)
  # calculate total length of each line
  clipped_layer$total_length <- st_length(clipped_layer$geometry)
  # convert from spatial polygon dataframe to dataframe without geometry
  shapefile_df <- as.data.frame(clipped_layer)
  # group by site and calculate density (and convert units to m/m^2)
  shapefile_density <- shapefile_df %>% 
    group_by(site) %>% 
    mutate(shapefile_density = sum(total_length)/catchment_area*10000) %>% 
    dplyr::select(site, catchment_area, shapefile_density) %>% 
    distinct()
  
  return(shapefile_density)
}

# function for human footprint index calculations
human_footprint_calcs <- function(human_footprint_data, spatial_extent_file, spatial_scale){
  # clip to intersection of spatial extent shapefile
  human_footprint_clip <- st_intersection(human_footprint_data, spatial_extent_file)
  # calculate area of each polygon
  human_footprint_clip$footprint_area <- st_area(human_footprint_clip$geometry)
  # convert from spatial polygon dataframe to dataframe without geometry
  human_footprint_df <- as.data.frame(human_footprint_clip)
  # group by site and impact intensity, and calculate percent area covered
  h_footrpint_percent <- human_footprint_df %>% 
    group_by(site, DN) %>% 
    mutate(percent_footprint = sum(footprint_area)/catchment_area*100) %>% 
    dplyr::select(site, catchment_area, DN, percent_footprint) %>% 
    rename(intensity = DN) %>% 
    distinct()
  # reshape into "wide" data and rename columns
  # removed 5 and added 3 (switch back if using 538m site buffer)
  if (spatial_scale == "small"){
    h_footprint_wide <- dcast(h_footrpint_percent, site~intensity, 
                              value.var = "percent_footprint") %>% 
      rename(human_impact_0 = "0", 
             human_impact_1 = "1", 
             human_impact_2 = "2", 
             human_impact_3 = "3", 
             human_impact_4 = "4", 
             # human_impact_5 = "5", 
             human_impact_6 = "6", 
             human_impact_7 = "7", 
             human_impact_8 = "8", 
             human_impact_9 = "9", 
             human_impact_10 = "10")
  }
  else{
    h_footprint_wide <- dcast(h_footrpint_percent, site~intensity, 
                              value.var = "percent_footprint", 
                              fun.aggregate=sum) %>% 
      rename(human_impact_0 = "0", 
             human_impact_1 = "1", 
             human_impact_2 = "2", 
             human_impact_3 = "3", 
             human_impact_4 = "4", 
             human_impact_5 = "5", 
             human_impact_6 = "6", 
             human_impact_7 = "7", 
             human_impact_8 = "8", 
             human_impact_9 = "9", 
             human_impact_10 = "10")
  }
  return(h_footprint_wide)
}

# function for disturbance calculations (repeat for each spatial scale)
disturbance_calculations <- function(spatial_extent_shapefile, spatial_scale){
  
  # calculate percent area for each disturbance and landscape layer
  insect_percent <- percent_of_catchment(insect, spatial_extent_shapefile) %>% 
    rename(insect_percent = percent_area)
  logging_percent <- percent_of_catchment(logging, spatial_extent_shapefile) %>% 
    rename(logging_percent = percent_area)
  cleared_percent <- percent_of_catchment(cleared, spatial_extent_shapefile) %>% 
    rename(cleared_percent = percent_area)
  fire_percent <- percent_of_catchment(fire, spatial_extent_shapefile) %>% 
    rename(fire_percent = percent_area)
  
  wetland_percent <- percent_of_catchment(wetland, spatial_extent_shapefile) %>% 
    rename(percent_wetland = percent_area)
  lakes_percent <- percent_of_catchment(lakes, spatial_extent_shapefile) %>% 
    rename(percent_lake = percent_area)
  barrens_percent <- percent_of_catchment(barrens, spatial_extent_shapefile) %>% 
    rename(percent_barren = percent_area)
  
  # calculate density of each layer
  paved_road_density <- density_in_catchment(paved_roads, spatial_extent_shapefile) %>% 
    rename(paved_road_density = shapefile_density)
  unpaved_road_density <- density_in_catchment(unpaved_roads, spatial_extent_shapefile) %>% 
    rename(unpaved_road_density = shapefile_density)
  trail_density <- density_in_catchment(trails, spatial_extent_shapefile) %>% 
    rename(trail_density = shapefile_density)
  
  # calculate % of each class of human footprint index
  human_footprint_wide <- human_footprint_calcs(human_footprint, 
                                                spatial_extent_shapefile, 
                                                spatial_scale)
  
  # combine into one dataframe
  if (spatial_scale == "small"){
    catchment_data <- site_list %>% 
      full_join((spatial_extent_shapefile %>% 
                   dplyr::select(site, catchment_area)), by="site") %>% 
      full_join(insect_percent, by=c("site", "catchment_area")) %>% 
      full_join(logging_percent, by=c("site", "catchment_area")) %>% 
      full_join(cleared_percent, by=c("site", "catchment_area"))  %>% 
      full_join(fire_percent, by=c("site", "catchment_area")) %>% 
      full_join(paved_road_density, by=c("site", "catchment_area")) %>% 
      full_join(unpaved_road_density, by=c("site", "catchment_area")) %>% 
      full_join(trail_density, by=c("site", "catchment_area")) %>% 
      full_join(barrens_percent, by=c("site", "catchment_area")) %>% 
      full_join(lakes_percent, by=c("site", "catchment_area")) %>% 
      full_join(wetland_percent, by=c("site", "catchment_area")) %>% 
      full_join(human_footprint_wide, by="site") %>% 
      mutate_at(vars(-park, -site, -geometry), as.numeric) %>% 
      replace(is.na(.), 0) %>% 
      mutate(total_disturbance = insect_percent + logging_percent + cleared_percent + fire_percent) %>% 
      mutate(total_road_density = paved_road_density + unpaved_road_density + trail_density) %>% 
      mutate(low_human_impact = human_impact_0) %>% 
      mutate(med_human_impact = human_impact_4 + human_impact_6)  %>% 
      mutate(high_human_impact = human_impact_8 + human_impact_10)
  }
  else{
    catchment_data <- site_list %>% 
      full_join((spatial_extent_shapefile %>% 
                   dplyr::select(site, catchment_area)), by="site") %>% 
      full_join(insect_percent, by=c("site", "catchment_area")) %>% 
      full_join(logging_percent, by=c("site", "catchment_area")) %>% 
      full_join(cleared_percent, by=c("site", "catchment_area"))  %>% 
      full_join(fire_percent, by=c("site", "catchment_area")) %>% 
      full_join(paved_road_density, by=c("site", "catchment_area")) %>% 
      full_join(unpaved_road_density, by=c("site", "catchment_area")) %>% 
      full_join(trail_density, by=c("site", "catchment_area")) %>% 
      full_join(barrens_percent, by=c("site", "catchment_area")) %>% 
      full_join(lakes_percent, by=c("site", "catchment_area")) %>% 
      full_join(wetland_percent, by=c("site", "catchment_area")) %>% 
      full_join(human_footprint_wide, by="site") %>% 
      mutate_at(vars(-park, -site, -geometry), as.numeric) %>% 
      replace(is.na(.), 0) %>% 
      mutate(total_disturbance = insect_percent + logging_percent + cleared_percent + fire_percent) %>% 
      mutate(total_road_density = paved_road_density + unpaved_road_density + trail_density) %>% 
      mutate(low_human_impact = human_impact_0 + human_impact_1 + human_impact_2 + human_impact_3) %>% 
      mutate(med_human_impact = human_impact_4 + human_impact_5 + human_impact_6)  %>% 
      mutate(high_human_impact = human_impact_7 + human_impact_8 + human_impact_9 + human_impact_10)
  }
  
  return(catchment_data) 
}

# use disturbance_calculations function to calculate data for each spatial scale
catchment_data <- disturbance_calculations(catchments, "large")
riparian_extent_data <- disturbance_calculations(riparian_extent, "medium")
site_extent_data <- disturbance_calculations(site_extent, "small")

# export data for each spatial scale as a csv
write.csv((catchment_data %>% dplyr::select(-geometry)), "C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc/output/disturbance_data_catchment.csv", row.names=FALSE)
write.csv((riparian_extent_data %>% dplyr::select(-geometry)), "C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc/output/disturbance_data_riparian.csv", row.names=FALSE)
write.csv((site_extent_data %>% dplyr::select(-geometry)), "C:/Users/hanna/OneDrive/Documents/GitHub/misc_msc/output/disturbance_data_local.csv", row.names=FALSE)

# 3) code for generating the "site" spatial extent (closest 10% of catchment area) ----
catchments_df <- catchments %>% 
  mutate(buffer_radius = sqrt((catchment_area*0.10)/pi)) %>%
  st_set_geometry(NULL)

# read in site pour points and generate buffer around each
# iterate through the data frame, adding points
pour_point <- st_read("site_pour_point.shp") %>% arrange(site)
stream_reach_df <- stream_reach %>% 
  st_set_geometry(NULL) %>% 
  mutate(buffer_radius = length*0.1) %>% 
  arrange(site)

buffer_sf <- st_buffer(pour_point$geometry, dist = as.numeric(catchments_df$buffer_radius), endCapStyle = "ROUND")
clipped_site_extent <- st_intersection(buffer_sf, catchments) 
proportional_site_extent <- st_intersection(clipped_site_extent, riparian_extent)
proportional_site_extent$catchment_area <- st_area(proportional_site_extent$geometry)

st_write(proportional_site_extent,"local_extent.shp", append=FALSE)
