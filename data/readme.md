# Empirical data files
*In situ* data collected from stream sites and geospatial data generated from satellite imagery in QGIS, used to find relationships between stream characteristics and terrestrial disturbances across spatial extents.  

## Folder directory  

### benthic_invertebrates.csv: 
Counts of benthic invertebrates collected from each stream site using a Surber sampler, identified to the family level by Entomogen Inc.  
|       Variable       |       Units          |      Description                                                                       |
|----------------------|----------------------|----------------------------------------------------------------------------------------|
| park                 | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova  |
| site                 | site name            | Unique code assigned to each stream site                                               |
| order                | taxonomic name       | Order of the benthic invertebrates counted                                             |
| family               | taxonomic name       | Family of the benthic invertebrates counted                                            | 
| count                | integer              | Count of benthic invertebrates in each taxonomic family                                |  

### canopy_cover.csv: 
Percent canopy cover each stream site, measured at 5 m intervals along the stream reach.  
|       Variable       |       Units          |      Description                                                                       |
|----------------------|----------------------|----------------------------------------------------------------------------------------|
| park                 | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova  |
| site                 | site name            | Unique code assigned to each stream site                                               |
| distance_upstream    | m                    | Distance upstream along the sampling reach that each canopy cover measurement was taken|
| canopy_cover         | %                    | Percent canopy cover along the cross section of the stream at each 5 m interval        |  

### channel_measurements.csv: 
Depth, width, and flow at each stream site, measured at three cross sections along the stream reach.
|       Variable       |       Units          |      Description                                                                                                               |
|----------------------|----------------------|--------------------------------------------------------------------------------------------------------------------------------|
| park                 | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova                                          |
| site                 | site name            | Unique code assigned to each stream site                                                                                       |
| cross_section        | integer              | Number assigned to the each of the three randomly selected cross seactions where channel measurements were taken at each site  |
| interval_across      | integer              | Interval across the stream width that depth and flow measurements were taken (number of measurements based on stream width)    |  
| depth                | m                    | Water depth at each interval across the stream channel                                                                         |
| flow1                | s                    | Elapsed time for the first flow measurement (used to calculate flow rate in m/s)                                               |  
| flow2                | s                    | Elapsed time for the second flow measurement (used to calculate flow rate in m/s)                                              | 
| flow3                | s                    | Elapsed time for the third flow measurement (used to calculate flow rate in m/s)                                               | 
| flow_distance        | m                    | Distance used for each timed flow measurement (used to calculate flow rate in m/s)                                             |
| wetted_width         | m                    | Width of the water in the stream channel                                                                                       |  
| bankful_width        | m                    | Width of the stream channel (between banks)                                                                                    | 

### chlorophyll_a.csv:
Spectrophotometer data from periphyton samples collected at each stream site. Absorbance values at key wavelengths are used to estimate periphyton biomass.
|       Variable                |       Units          |      Description                                                                        |
|-------------------------------|----------------------|-----------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova   |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such     |
| sample_number                 | integer              | Number assigned to individual subsamples of the combined sample from each stream site   |
| absorbance_664nm_unacidified  | positive value       | Value of the absorbance measurement taken at the 664 nm wavelength before acidification |
| absorbance_750nm_unacidified  | positive value       | Value of the absorbance measurement taken at the 750 nm wavelength before acidification |
| absorbance_664nm_acidified    | positive value       | Value of the absorbance measurement taken at the 664 nm wavelength after  acidification |
| absorbance_750nm_acidified    | positive value       | Value of the absorbance measurement taken at the 750 nm wavelength after acidification  |

### doc.csv:  
Dissolved organic carbon data from filtered water samples collected at each stream site, measured with a DOC/TDN analyzer.
|       Variable                |       Units          |      Description                                                                      |
|-------------------------------|----------------------|---------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such   |
| sample_number                 | integer              | Number assigned to individual subsamples of the combined sample from each stream site |
| date                          | text                 | Date and time the absorbance measurements were taken                                  |
| injection_number              | integer              | Number assigned to individual measurements taken from each subsample                  |
| area                          | positive number      | Area under the curve of the plot generated by the TOC analyzer, used to calculate DOC |
| doc                           | mg/L                 | Value of the DOC measurement                                                          |

### tn.csv:  
Dissolved nitrogen data from filtered water samples collected at each stream site, measured with a DOC/TDN analyzer.  
|       Variable                |       Units          |      Description                                                                      |
|-------------------------------|----------------------|---------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such   |
| sample_number                 | integer              | Number assigned to individual subsamples of the combined sample from each stream site |
| date                          | text                 | Date and time the absorbance measurements were taken                                  |
| injection_number              | integer              | Number assigned to individual measurements taken from each subsample                  |
| area                          | positive number      | Area under the curve of the plot generated by the TOC analyzer, used to calculate TN  |
| tn                            | mg/L                 | Value of the TN measurement                                                           |

### pebble_count.csv:  
Counts of substrate size and embeddedness at each stream site, following Canadian Aquatic Biomonitoring Network (CABIN) guidelines (CABIN Field Manual, 2009).
|       Variable                |       Units          |      Description                                                                       |
|-------------------------------|----------------------|----------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova  |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such    |
| pebble_number                 | integer              | Number assigned to individual pebbles measured at each stream site                     |
| intermediate_axis             | positive value       | Width of the intermediate axis of each pebble                                          |
| embeddedness                  | %                    | Percent embeddedness of every 10th pebble                                             |

### periphyton_foil.csv:  
Mass of foil used to cover the surface area of the rocks that the periphyton samples were collected from. These values were converted to surface area following (Hauer & Lamberti, 2007).
|       Variable                |       Units          |      Description                                                                                    |
|-------------------------------|----------------------|-----------------------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova               |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such                 |
| date                          | text                 | Date the periphyton samples were collected                                                          |
| sample_number                 | integer              | Number assigned to individual periphyton samples                                                    |
| foil_mass                     | g                    | Mass of the foil required to covre each rock samples for periphyton, used to calculate area sampled |

### surber_sampling.csv:  
Number of surber samples collected from each stream site.  
|       Variable                |       Units          |      Description                                                                                    |
|-------------------------------|----------------------|-----------------------------------------------------------------------------------------------------|
| park                          | GM or TN             | Provincial park the site is within or closest to; GM = Gros Morne and TN = Terra Nova               |
| site                          | site name            | Unique code assigned to each stream site; standards and blanks are labelled as such                 |
| surber_samples                | integer              | Number of Surber samples collectected at each stream site, used to measure taxa per area            |

### tss_filters.csv:  
Mass of total suspended solids measured from water samples collected at each stream site. Note that these measurements were below the instrument detection limit and were not included in our analyses.  

### water_chemistry.csv:  
Measurements of ph, water temperature, electrical conductivity, total dissolved solids, alkalinity, and turbidity from each stream site.  

### periphyton_afdm.csv:  
Ash free dry mass (AFDM) of periphyton samples collected at each stream site.  

### spatial_data.zip:  
Zip file containing shapefiles used to calculate disturbance metrics at each site.  
* stream_reach.shp: line showing the sampling reach at each stream site  
* stream_sites.shp: points of sampling location for each stream site
* catchments.shp: polygons of the catchments upstream of the sampling location at each stream site (largest spatial extent)
* riparian_extent.shp: polygons of the mid-sized spatial extent (100 m riparian buffer on each side of the stream and tributaries)  
* local_extent.shp: polygons of the smallest spatial extent at each stream site (10% of the catchment area, closest to the sampling location)  
* forest_disturbance.shp: polygons of forest disturbance from logging, insect outbreaks, fores fire, and a general "cleared" category within the site catchments  
* paved_roads.shp: lines of all paved roads within the site catchments
* trails.shp: lines of all trails within the site catchments  
* unpaved_roads.shp: lines of all unpaved roads (including ATV trails) within the site catchments  
