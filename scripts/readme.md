# Code for processing empirical data and simulating disturbances in our meta-ecosystem model  
This folder contains the main scripts for 1) processing and 2) analyzing our empirical data, and 3) simulating disturbance in our meta-ecosystem model.

1. Stream data os cleaned using the *stream_data_cleaning.R* script, and disturbance metrics are calculated at three spatial extents at each stream site using *calculating_disturbance_metrics.R*
2. The *correlogram.R* script is used to make sure that there is no significant correlation between variables used in our statistical models, and the *stats_models.R* script is used to compare general linear models for each of our stream quality metrics to determine relationships between terrestrial disturbance and stream quality.
3. The *meta_ecosystem_analytical_equilibria.nb* script is used to generate analytical equilibria and a jacobian matrix for the meta-ecosystem model, which is them used in the *simulate_disturbance.R* script to first generate 10,000 stable equilibria and simulate two types of terrestrial disturbance across a range of intensities. The *global_sensitivity_analysis.R* script was used during this process determine key parameters that could be used to simulate disturbance in our model.

Sections of code requiring data files that are not provided are commented out, but can be added back in as needed.  

## Scripts  
### Processing empirical data  
* stream_data_cleaning.R  
    requires: water_chemistry.csv, canopy_cover.csv, pebble_count.csv, channel_measurements.csv, chlorophyll_a.csv, periphyton_foil.csv, periphyton_afdm.csv
    outputs:  
* calculating_disturbance_metrics.R  
    requires:  
    outputs:  

### Analyzing empirical data  
* correlogram.R
    requires:  
    outputs: 
* stats_models.R
    requires:  
    outputs: 

### Meta-ecosystem model  
* global_sensitivity_analysis.R
    requires:  
    outputs: 
* meta_ecosystem_analytical_equilibria.nb
    requires:  
    outputs: 
* simulate_disturbance.R
    requires:  
    outputs: 

## Software and packages

## References  
