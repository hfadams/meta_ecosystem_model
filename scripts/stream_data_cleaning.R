# formatting stream data and calculating new metrics:
# 1) Water chemistry
# 2) Channel morphology
# 3) Canopy cover
# 4) Periphyton
# 5) Benthic invertebrates
# 6) Compile all data together

## load packages ----
library(dplyr)
library(ggplot2)
library(lattice)
library(stringr)
library(reshape2)
library(asbio)
library(vegan)

## read in files ----
water_chem <- read.csv("data/water_chemistry.csv")  %>% 
  rename("conductivity" = ec)
canopy <- read.csv("data/canopy_cover.csv")
pebble_count <- read.csv("data/pebble_count.csv")
channel <- read.csv("data/channel_measurements.csv")
chlorophyll <- read.csv("data/chlorophyll_a.csv")
periphyton_foil <- read.csv("data/periphyton_foil.csv")
periphyton_afdm <- read.csv("data/periphyton_afdm.csv")
tss <- read.csv("data/tss_filters.csv")
inverts <- read.csv("data/benthic_invertebrates.csv")
doc <- read.csv("data/doc.csv")
tn <- read.csv("data/tn.csv")
surber_sample_area <- read.csv("data/surber_sampling.csv")
invert_traits <- read.csv("data/invert_traits_usgs.csv") %>% 
  dplyr::select(Family, 
                Measured_length, 
                Feed_mode_prim, 
                Feed_mode_sec, 
                Feed_mode_comments) %>% 
  rename("family" = Family)
invert_traits$Measured_length <- as.numeric(invert_traits$Measured_length)
invert_coeficients <- read.csv("data/invert_coeficients.csv")

## calculate site metrics:

# 1) water chemistry ----

# calculate mean and sd for water samples
water_chem_means <- water_chem %>% 
  select(park,
         site,
         ph,
         temperature,
         conductivity,
         tds,
         alkalinity) %>%
  group_by(site) %>%
  summarise(across(ph:alkalinity, ~mean(.x, na.rm = TRUE)))

water_chem_sd <- water_chem %>% 
  select(park,
         site,
         ph,
         temperature,
         conductivity,
         tds,
         alkalinity) %>%
  group_by(site) %>%
  summarise(across(ph:alkalinity, ~sd(.x, na.rm = TRUE)))

# total suspended solids filters --
# first, calculate limit of detection (LOD) from filter blanks
tssblk1 = 0.104
tssblk2 = 0.103
tssblk3 = 0.102

tss_lod = 2*sd(c(tssblk1, tssblk2, tssblk3))
mean_filter_mass = mean(c(tssblk1, tssblk2, tssblk3))

# filter for samples below the LOD and calculate mean (some samples required multiple filters)
tss_calcs <- tss %>% 
  filter(dry_mass > (mean_filter_mass+tss_lod)) %>%
  mutate(tss_mg_l = (dry_mass - mean_filter_mass)*100/sample_volume*1000) %>% 
  group_by(site) %>% 
  mutate(tss = mean(tss_mg_l)) %>% 
  select(site, 
         tss) %>% 
  distinct()
# not enough samples above limit of detection, cannot use tss data!

# DOC and TDN --
# calculate method detection limit (MDL) for doc and tdn based on field blanks
mdl_doc <- doc %>% 
  filter(sample_number == "blank") %>% 
  mutate(doc = 2*sd(doc_mgl)) %>% 
  select(doc) %>% 
  distinct()
samples_under_doc_mdl <- doc %>% 
  filter(doc_mgl < mdl_doc$doc) %>% 
  filter(sample_number != "blank")

mdl_tn <- tn %>% 
  filter(sample_number == "blank") %>% 
  mutate(tn = 2*sd(tn_mgl)) %>% 
  select(tn) %>% 
  distinct()
samples_under_tn_mdl <- tn %>% 
  filter(tn_mgl < mdl_tn$tn) %>% 
  filter(sample_number != "blank")

# filter out rows where samples failed QA/QC (GM-MAS-3 tdn and doc), GM LOM 1 tdn injection # 1 and 2
doc_calcs <- doc %>% 
  filter(sample_number != "blank") %>%
  slice(-c(226,227,228)) %>% 
  group_by(site) %>% 
  mutate(doc = mean(doc_mgl), doc_sd = sd(doc_mgl)) %>% 
  ungroup() %>% 
  select(site, 
         doc, 
         doc_sd) %>% 
  distinct()

tn_calcs <- tn %>% 
  filter(sample_number != "blank") %>%
  slice(-c(240,241,242, 204, 205)) %>% 
  group_by(site) %>% 
  mutate(tn = mean(tn_mgl), tn_sd = sd(tn_mgl)) %>% 
  ungroup() %>% 
  select(site, tn, tn_sd) %>% 
  distinct()
# GM-CCR-2 is flagged as higher than expected, but leaving in the dataset)


# 2) Channel morphology and stream flow ----

# i) pebble count
# Make a function to assign class based on substrate size
rock_class <- function(dataframe, axis_length){
  class <- case_when(axis_length < 0.1 ~ "fine_sand_silt_clay",
                     axis_length < 0.2 ~ "course_sand",
                     axis_length < 1.6 ~ "gravel",
                     axis_length < 3.2 ~ "small_pebble",
                     axis_length < 6.4 ~ "large_pebble",
                     axis_length < 12.8 ~ "small_cobble",
                     axis_length < 25.6 ~ "large_cobble",
                     axis_length < 100 ~ "boulder")
  
  # append to the original data frame
  dataframe$pebble_class = class
  return(dataframe)
}

# make a function to do all the calculations on pebble_count
pebble_calculations <- function(data){
  # Calculate the Wolman D50 (mean diameter) and assign class based on this mean
  pebble_means <- data %>% 
    group_by(site) %>% 
    summarise(across(intermediate_axis:embeddedness, ~mean(.x, na.rm = TRUE))) %>%
    rename(wolmanD50=intermediate_axis)
  
  pebble_sd <- data %>% 
    group_by(site) %>% 
    summarise(across(intermediate_axis:embeddedness, ~sd(.x, na.rm = TRUE))) %>%
    rename(wolmanD50_sd=intermediate_axis, embeddedness_sd=embeddedness)
  
  classified_pebbles <- rock_class(data, data$intermediate_axis)
  classified_pebbles_max <- classified_pebbles %>% 
    group_by(site) %>% 
    mutate(num_pebbles=max(pebble_number))
  
  percent_comp <- classified_pebbles_max %>% 
    group_by(site, pebble_class) %>% 
    mutate(class_count=length(pebble_class)) %>% 
    mutate(class_percent=(class_count/num_pebbles*100)) %>% 
    select(park, 
           site, 
           pebble_class, 
           class_count, 
           class_percent) %>% 
    distinct()
  
  # compile into one dataframe based on site
  all_pebble_data <- percent_comp %>% 
    merge(pebble_means, by="site") %>% 
    merge(pebble_sd, by="site")
  
  return(all_pebble_data) 
}

all_pebble_data <- pebble_calculations(pebble_count)

# select columns to later export
pebble_calcs <- all_pebble_data %>% 
  select(site, 
         wolmanD50, 
         wolmanD50_sd,
         embeddedness,
         embeddedness_sd) %>% 
  distinct()

# 1) flow
# calculate mean flow from the three replicates, then means for each site
channel_selected <- channel %>% 
  mutate(flow = (flow1 + flow2 + flow3)/flow_distance/3) %>% 
  select(-flow1,
         -flow2,
         -flow3,
         -flow_distance,
         -notes)

channel_means <- channel_selected %>%
  group_by(site) %>%
  summarise(across(depth:flow, ~mean(.x, na.rm = TRUE)))

channel_sd <- channel_selected %>%
  group_by(site) %>%
  summarise(across(depth:flow, ~sd(.x, na.rm = TRUE)))

# 3) canopy cover ----
# make columns for reach length and average cover
canopy_calcs <- canopy %>% 
  group_by(site) %>% 
  mutate(reach_length = max(distance_upstream), 
         canopy = mean(canopy_cover),
         canopy_sd = sd(canopy_cover)) %>% 
  select(park,
         site, 
         reach_length, 
         canopy,
         canopy_sd) %>% 
  distinct()

# 4) Periphyton ----
# first, calculate the area of the rocks sampled
foil_mass_cm2 <- (0.614 + 0.656)/2/10  # mean mass of 1 cm^2 of aluminum foil
periphyton_area <- periphyton_foil %>% 
  mutate(rock_area_cm2 = foil_mass/foil_mass_cm2)

# a) calculate periphyton biomass using ash free dry mass method, and filter --
# for samples over the limit of detection (same as tss since using same balance)
peri_lod = tss_lod

peri_afdm <- periphyton_afdm %>% 
  full_join(periphyton_area, by=c("park", 
                                  "site", 
                                  "sample_number")) %>% 
  group_by(site) %>% 
  mutate(total_rock_area = sum(rock_area_cm2, na.rm = TRUE)) %>% 
  mutate(biomass_g = (dry_mass-afdm)) %>% 
  filter(biomass_g > peri_lod) %>%
  mutate(biomass_mg = biomass_g*100, total_biomass = sum(biomass_mg, na.rm = TRUE)) %>% 
  mutate(biomass_cm2 = total_biomass/total_rock_area*2) %>% # multiply by 2 because half filter
  select(date, 
         site, 
         total_rock_area, 
         total_biomass, 
         biomass_cm2) %>% 
  distinct()

# 18 filter samples were below LOD, will not use afdm data in the analysis


# b) calculate periphyton biomass using chlorophyll a --

path_length <- 1
extraction_volume <- 7 

# calculate method detection limit
chlorophyll_blanks = chlorophyll %>% 
  filter(site == "blank") %>%
  filter(sample_number != 4) %>% # mislabeled standard
  mutate(chlorophyll_concentration = 11*(absorbance_664nm_unacidified-absorbance_750nm_unacidified)*extraction_volume/path_length) %>%
  summarise(across(chlorophyll_concentration, ~sd(.x, na.rm = TRUE)))

chlorophyll_mdl <- 2*chlorophyll_blanks$chlorophyll_concentration

# calculate chlorophyll a per unit area (*2  because half filter)
# filter for samples above mdl

chlorophyll_prelim_calcs <- chlorophyll %>% 
  inner_join(periphyton_area, by=c("park", 
                                   "site", 
                                   "sample_number")) %>%
  mutate(chlorophyll_concentration = 11*(absorbance_664nm_unacidified-absorbance_750nm_unacidified)*extraction_volume/path_length) %>%
  filter(chlorophyll_concentration > chlorophyll_mdl) %>%
  mutate(periphyton_cm2 = chlorophyll_concentration/rock_area_cm2*2)

# calculate mean and standard deviation for each site
chlorophyll_calcs <- chlorophyll_prelim_calcs %>% 
  group_by(site) %>% 
  mutate(periphyton_biomass = mean(periphyton_cm2, na.rm = TRUE)) %>% 
  mutate(periphyton_biomass_sd = sd(periphyton_cm2, na.rm = TRUE)) %>% 
  select(site, 
         periphyton_biomass,
         periphyton_biomass_sd) %>% 
  distinct()


# 5) benthic invertebrates ----
surber_area <- 30.48**2 # area sampled in square cm (12"x12")
inverts[is.na(inverts)] = 0  # replace nan with 0
inverts_unique <- inverts %>% 
  select(order, 
         family) %>% 
  distinct()

# calculate number of individuals, then sift the dataframe to count the number 
# of unique taxa and EPT index
invert_calcs <- inverts %>% 
  group_by(site) %>% 
  mutate(total_individuals = sum(invert_count)) %>% 
  filter(invert_count != 0) %>% 
  group_by(site) %>% 
  mutate(unique_taxa = length(invert_count)) %>% 
  filter(order == "ephemeroptera" | order == 
           "plecoptera" | order == "trichoptera") %>% 
  mutate(ept_index = sum(invert_count)/total_individuals*100) %>% 
  select(park, 
         site, 
         total_individuals, 
         unique_taxa, 
         ept_index) %>% 
  distinct()

# merge with sampling area dataframe and calculate individuals and taxa per area
inverts_per_area <- surber_sample_area %>% 
  full_join(invert_calcs, by=c("site", 
                               "park")) %>% 
  mutate(invertebrates_cm2 = total_individuals/(surber_area*surber_samples), 
         taxa_cm2 = unique_taxa/(surber_area*surber_samples)) %>% 
  select(site, 
         surber_samples, 
         total_individuals, 
         unique_taxa, 
         ept_index, 
         invertebrates_cm2, 
         taxa_cm2)

# calculate invertebrate biomass
# extract average length for each order
invert_traits_mean <- left_join((inverts %>% 
                                   select(order, 
                                          family) %>% 
                                   distinct()), 
                                invert_traits, 
                                by="family") %>% 
  group_by(order) %>% 
  mutate(mean_length = mean(Measured_length, na.rm=TRUE)) %>% 
  select(order, 
         mean_length) %>% 
  distinct()

# add length for orders with missing data points
invert_traits_mean$mean_length[invert_traits_mean$order == "collembola"] <- 1.28
invert_traits_mean$mean_length[invert_traits_mean$order == "bivalvia"] <- 4.6
invert_traits_mean$mean_length[invert_traits_mean$order == "gastropoda"] <- 8.32

# determine functional groups for each family
functional_groups <- inverts %>% 
  select(order, 
         family) %>% 
  distinct() %>% 
  left_join(invert_traits, by="family") %>% 
  filter(!is.na(Feed_mode_prim)) %>% 
  filter(Feed_mode_prim != "") %>% 
  group_by(order, family) %>%  
  filter(row_number()==1) %>% 
  select(!Measured_length) %>% 
  right_join(inverts_unique, by=c("order", "family")) %>% 
  rename(functional_group = Feed_mode_prim)

# add functional group for families where data is missing
functional_groups$functional_group[functional_groups$family == "Collembola"] <- "Collector_filterer"
functional_groups$functional_group[functional_groups$family == "immature"] <- "Shredder" # most common group for most taxa
functional_groups$functional_group[functional_groups$family == "Blebphariceridae"] <- "Collector_gatherer"
functional_groups$functional_group[functional_groups$family == "Hirudinea"] <- "Predator"
functional_groups$functional_group[functional_groups$family == "Sarcoptiformes"] <- "Predator"
functional_groups$functional_group[functional_groups$family == "Noctuidae"] <- "Scraper_grazer" # terrestrial; shouldn't be counted?
functional_groups$functional_group[functional_groups$family == "Leptophlebiidae"] <- "Scraper_grazer"
functional_groups$functional_group[functional_groups$family == "Odontoceridae"] <- "Shredder"
functional_groups$functional_group[functional_groups$family == "Elmidae"] <- "Scraper_grazer"
functional_groups$functional_group[functional_groups$family == "Enchytraeidae"] <- "Collector_gatherer"
functional_groups$functional_group[functional_groups$family == "Lumbriculidae"] <- "Collector_gatherer"
functional_groups$functional_group[functional_groups$family == "Naididae"] <- "Collector_gatherer"
functional_groups$functional_group[functional_groups$family == "Hyalellidae"] <- "Collector_gatherer"
functional_groups$functional_group[functional_groups$family == "Sisyridae"] <- "Predator"
functional_groups$functional_group[functional_groups$family == "Gyrinidae"] <- "Scavenger"
functional_groups$functional_group[functional_groups$family == "Pisidiidae"] <- "Collector_filterer"

# replace small feeding group categories
functional_groups <- functional_groups %>% 
  mutate(across("functional_group", str_replace, "Scavenger", "collector_gatherer")) %>% 
  mutate(across("functional_group", str_replace, "Parasite", "predator")) %>% 
  mutate(across("functional_group", str_replace, "Collector-gatherer", "collector_gatherer")) %>% 
  mutate(across("functional_group", str_replace, "Collector-filterer", "collector_filterer")) %>% 
  mutate(across("functional_group", str_replace, "Scraper/grazer", "scraper_grazer")) %>% 
  mutate(across("functional_group", str_replace, "Shredder", "shredder")) %>% 
  mutate(across("functional_group", str_replace, "Collector_gatherer", "collector_gatherer")) %>% 
  mutate(across("functional_group", str_replace, "Collector_filterer", "collector_filterer")) %>% 
  mutate(across("functional_group", str_replace, "Scraper_grazer", "scraper_grazer")) %>% 
  mutate(across("functional_group", str_replace, "Predator", "predator")) %>% select(family,
                                                                                     order,
                                                                                     functional_group)

# Add "all insect" coefficients where data is not available
# Use allometric equation (power law) to estimate biomass based on length
invert_biomass_feeding <- invert_traits_mean %>% 
  left_join(invert_coeficients, by="order") %>% 
  mutate_at(vars(b), ~replace(., is.na(.), 2.788)) %>% 
  mutate_at(vars(a), ~replace(., is.na(.), 0.0064)) %>% 
  mutate(biomass_individual = a*mean_length**b) %>% 
  right_join(inverts, by="order") %>% 
  mutate(biomass_order = biomass_individual*invert_count) %>% 
  group_by(site) %>% 
  mutate(invert_biomass_mg=sum(biomass_order)) %>% 
  left_join(surber_sample_area, by=c("park", 
                                     "site")) %>% 
  mutate(invertebrate_biomass = invert_biomass_mg/(surber_area*surber_samples)) 

ept_biomass <- invert_biomass_feeding %>% 
  group_by(site) %>% 
  filter(order == "ephemeroptera" | order == "plecoptera" | order == "trichoptera") %>% 
  mutate(ept_biomass_cm2 = sum(biomass_order)/(surber_area*surber_samples)) %>% 
  select(site,
         ept_biomass_cm2) %>% 
  distinct()

predator_biomass <- invert_biomass_feeding %>% 
  full_join(functional_groups) %>% 
  group_by(site) %>%
  filter(functional_group == "predator") %>% 
  mutate(predator_biomass_cm2 = sum(biomass_order)/(surber_area*surber_samples)) %>% 
  select(site,
         predator_biomass_cm2) %>% 
  distinct()

invert_biomass_data <- invert_biomass_feeding %>% 
  full_join(ept_biomass, by="site") %>%
  full_join(predator_biomass, by="site") %>% 
  select(site,
         invertebrate_biomass,
         ept_biomass_cm2,
         predator_biomass_cm2) %>% 
  distinct()

# calculate counts of functional groups
invert_fun_groups <- inverts %>% 
  left_join(functional_groups, by=c("order", "family")) %>% 
  group_by(site, functional_group) %>% 
  mutate(fun_group_count=sum(invert_count)) %>% 
  select(site, functional_group, fun_group_count) %>% 
  distinct() %>% 
  left_join((invert_calcs %>% select(site, total_individuals)), by="site") %>% 
  mutate(percent_fun_group = fun_group_count/total_individuals*100)
invert_fun_groups_wide <- dcast(invert_fun_groups, site~functional_group, 
                                value.var = "percent_fun_group")

# 6) Compile data by site and export ----

# compile mean data
all_mean_data <- canopy_calcs %>%
  select(park, site, canopy, reach_length) %>%
  full_join(water_chem_means, by = "site") %>% 
  full_join((doc_calcs %>% select(site, doc)), by="site") %>%
  full_join((tn_calcs %>% select(site, tn)), by="site") %>%
  full_join((pebble_calcs %>% select(site, wolmanD50, embeddedness)), by="site") %>% 
  full_join(channel_means, by="site") %>% 
  full_join((chlorophyll_calcs %>% select(site, periphyton_biomass)), by="site") %>% 
  full_join((invert_calcs %>% select(site, ept_index)), by="site") %>%
  full_join((invert_biomass_data %>% select(site, invertebrate_biomass)), by="site") %>% 
  full_join(invert_fun_groups_wide, by="site")
  

# compile standard deviation data
all_sd_data <- canopy_calcs %>%
  select(park, site, canopy_sd) %>%
  full_join(water_chem_sd, by = "site") %>% 
  full_join((doc_calcs %>% select(site, doc_sd)), by="site") %>%
  full_join((tn_calcs %>% select(site, tn_sd)), by="site") %>%
  full_join((pebble_calcs %>% select(site, wolmanD50_sd, embeddedness_sd)), by="site") %>% 
  full_join(channel_sd, by="site") %>% 
  full_join((chlorophyll_calcs %>% select(site, periphyton_biomass_sd)), by="site") %>%
  rename("doc"=doc_sd, 
         "tn"=tn_sd,
         "canopy"=canopy_sd,
         "wolmanD50"=wolmanD50_sd,
         "embeddedness"=embeddedness_sd,
         "periphyton_biomass"=periphyton_biomass_sd)
  
  
# export files
write.csv(all_mean_data, "output/empirical_stream_data.csv", row.names = FALSE)
write.csv(all_sd_data, "output/empirical_stream_data_standard_deviations.csv", row.names = FALSE)
