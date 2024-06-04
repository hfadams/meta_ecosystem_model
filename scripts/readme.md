# Code for processing empirical data and simulating disturbances in the meta-ecosystem model  
This folder contains the main scripts for processing and analyzing our empirical data, and simulating disturbance in our meta-ecosystem model:

1. Stream data is cleaned using the *stream_data_cleaning.R* script, and disturbance metrics are calculated at three spatial extents at each stream site using *calculating_disturbance_metrics.R*
2. The *correlogram.R* script is used to make sure that there is no significant correlation between variables used in our statistical models, the *stats_models.R* script is used to compare general linear models for each of our stream quality metrics and calculate AICc, and *model_averaging.r* is used to perform model aaveraging after removing models that have uninformative parameters and are within delta AICc 2 of the null model. The results are used to determine relationships between terrestrial disturbance and stream quality.
3. The *meta_ecosystem_analytical_equilibria.nb* script is used to generate analytical equilibria and a jacobian matrix for the meta-ecosystem model, which is them used in the *simulate_disturbance.R* script to first generate 10,000 stable equilibria and simulate two types of terrestrial disturbance across a range of intensities. The *global_sensitivity_analysis.R* script was used during this process determine key parameters that could be used to simulate disturbance in our model.

Sections of code requiring data files that are not provided are commented out, but can be added back in as needed.  

## Scripts  
### Processing empirical data  
* stream_data_cleaning.R  
    requires: water_chemistry.csv, canopy_cover.csv, pebble_count.csv, channel_measurements.csv, chlorophyll_a.csv, periphyton_foil.csv, periphyton_afdm.csv
    outputs: all_mean_data.csv, all_sd_data.csv (mean and standard deviations for *in situ* data from each stream site)  
* calculating_disturbance_metrics.R  
    requires: stream_reach.shp, catchments.shp, stream_sites.shp, riparian_extent.shp, local_extent.shp, paved_roads.shp, unpaved_roads.shp, trails.shp, forest_disturbance.shp, lakes.shp, barrens.shp, wetland.shp, human_footrpint.shp  
    outputs: disturbance_data_local.csv, disturbance_data_riparian.csv, disturbance_data_catchment.csv  

### Analyzing empirical data  
* correlogram.R  
    requires: disturbance_data_catchment.csv  
    outputs: disturbance_correlogram.svg  
* stats_models.r  
    requires: disturbance_data_local.csv, disturbance_data_riparian.csv, disturbance_data_catchment.csv  
    outputs: empirical_glm_results.csv, rescaled_data_local.csv, rescaled_data_riparian.csv, rescaled_data_catchment.csv    
* model_averaging.r
    requires: rescaled_data_local.csv, rescaled_data_riparian.csv, rescaled_data_catchment.csv (generated from stats_models.r)  
    outputs: model_summary_local_averaged.svg, model_summary_local_averaged.svg, model_summary_local_averaged.svg, model_averaging_formatted.csv  

### Meta-ecosystem model  
* meta_ecosystem_analytical_equilibria.nb
    requires: NA   
    outputs: analytical equilibria and jacobian matrix
* simulate_disturbance.R
    requires: NA 
    outputs: param_estimates1.csv, param_estimates2.csv, param_simulations_stable.csv, stable_data_θt_μt.csv, stable_data_βa_αa.csv, surface_data_forestry.csv, surface_data_atv.csv, surafce_plots_forestry.svg, surafce_plots_nutrients_forestry.svg, surafce_plots_atv.svg, surafce_plots_atv_nutrients.svg  
* global_sensitivity_analysis.R
    requires: param_simulations_stable_10000.csv  
    outputs: gsa_plot.svg, gsa_results.csv  

## Packages
library(dplyr)
library(ggplot2)
library(lattice)
library(stringr)
library(reshape2)
library(asbio)
library(vegan)
library(raster)
library(rgdal)
library(sf)
library(elevatr)
library(lwgeom)
library(GGally)
library(ggplot2)
library(svglite)
library(gridExtra)
library(tibble)
library(wiqid)
library(nlme)
library(AICcmodavg)
library(deSolve)
library(rootSolve)
library(tidyverse)
library(ggpol) 
library(MetBrewer) 
library(ggpubr)
library(randomForest)
library(lhs)
library(data.table)

## References  
