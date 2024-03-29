# Meta-ecosystem model

The scripts within this repository are used to 1) clean and analyze empirical data collected from streams in Gros Morne National Park and Terra Nova National Park in Newfoundland, Canada and 2) simulate terrestrial disturbance in a meta-ecosystem model of these stream-riparian ecosystems. These files are inteded to make our analyses as transparent and reproducible as possible, and have the stream data openly accessible for further analysis. A thorough description of the collection methods can be found in the appendix of the associated manuscript entitled "Integrating field data and a meta-ecosystem model to study the effects of multiple terrestrial disturbances on small stream ecosystem function". 

The land class data are not publically available, but Forest Resource Inventory shapefiles can be requested from the Newfoundland government, or can be digitized from satellite imagery.

## Data sources  

1) Field collection July-August 2022   
Water quality and stream measurements, found in "data" folder  

2) United States Geological Survey (USGS) invertebrate database  
Vieira, N.K.M. et al. (2016) ‘A Database of Lotic Invertebrate Traits for North America’, U.S. Geological Survey Data Series, 187, pp. 1–15.  
[USGS Database](https://doi.org/10.3133/ds187)  

    Read their data sharing policy [here](https://www.usgs.gov/media/files/casc-data-sharing-policy).

3) Benthic invertebrate mass to length conversion coefficients
   Benke, A. C., Huryn, A. D., Smock, L. A., & Wallace, J. B. (1999). Length-Mass Relationships for Freshwater Macroinvertebrates in North America with Particular Reference to the Southeastern United States. In Source: Journal of the North American Benthological Society (Vol. 18, Issue 3).
   
4) CanElevation 5m resolution digital elevation model  
Government of Canada (2022) High Resolution Digital Elevation Model (HRDEM) - CanElevation Series.
[CanElevation HRDEM](https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995)

5) Human impact index shapefile  
   Government of Newfoundland and Labrador. (2013). Newfoundland and Labrador Human Footprint: A Snapshot of Human Influence on the Landscape. Department of Environment and Conservation, 25.  

6) Data sharing agreement with Gros Morne National Park and Terra Nova National Park.
    
7) Data sharing agreement with the Government of Newfoundland and Labrador.

## Methods  
### In situ data  
We collected the in situ data in triplicate at each stream site and took the mean of the three samples. Refer to the manuscript and appendix A for detailed methods.

### Model simulations
After developing a riparian-stream meta-ecosystem model, we used Mathematica to solve for all possible analytical equilibria of the model. We then selected the equilibrium that was locally stable and feasible, and parameterized the model by generating 10,000 random parameter combinations, each within a range from 0-10 (or 0-1 if a proportion). From these we selected the first 1,000 equilibria taht were feasible, locally stable, and where the benthic invertebrate biomass was greater than periphyton biomass. We used these simulations as the "undisturbed" meta-ecosystem to which we created "terrestrial disturbances" by increasing key parameters to simulate tree removal and increased erosion. Refer to the manuscript and Appendix B for further details.

### Quality assurance  
We removed all data below the method detection limit of each *in situ* measurement and lab analysis before statistical analysis. We ensured that there was no correlation between predictor variables in the empirical dataset before developing the general linear models. 

We performed a global sensitivity analysis on each trophic level and productivity metric in the meta-ecosystem model to identify parameters creating the most uncertainty in the model.  

## software and packages
All data processing and analyses for this project were implemented using R (ver. 4.2.2), Mathematica (ver. 13.2.1), and QGIS (ver. 3.26.3).

## repository directory
### Folder 1: data
Empirical data used for statistical analysis (*in situ* data collected from stream sites and shapefiles digitized in QGIS)  
* [benthic_invertebrates.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/benthic_invertebrates.csv): Counts of benthic invertebrates collected from each stream site using a Surber sampler, identified to the family level by Entomogen Inc.
* [canopy_cover.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/canopy_cover.csv): Percent canopy cover each stream site, measured at 5 m intervals along the stream reach.
* [channel_measurements.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/channel_measurements.csv): Depth, width, and flow at each stream site, measured at three cross sections along the stream reach.
* [chlorophyll_a.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/chlorophyll_a.csv): Spectrophotometer data from periphyton samples collected at each stream site. Absorbance values at key wavelengths are used to estimate periphyton biomass.
* [doc.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/doc.csv): Dissolved organic carbon data from filtered water samples collected at each stream site, measured with a DOC/TDN analyzer.
* [pebble_count.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/pebble_count.csv): Counts of substrate size and embeddedness at each stream site, following Canadian Aquatic Biomonitoring Network (CABIN) guidelines (CABIN Field Manual, 2009).
* [periphyton_foil.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/periphyton_foil.csv): Mass of foil used to cover the surface area of the rocks that the periphyton samples were collected from. These values were converted to surface area following (Hauer & Lamberti, 2007).
* [spatial_data.zip](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/spatial_data.zip): zip file containing shapefiles used to calculate disturbance metrics at each site
    * stream_reach.shp: line showing the sampling reach at each stream site  
    * stream_sites.shp: points of sampling location for each stream site
    * catchments.shp: polygons of the catchments upstream of the sampling location at each stream site (largest spatial extent)
    * riparian_extent.shp: polygons of the mid-sized spatial extent (100 m riparian buffer on each side of the stream and tributaries)  
    * local_extent.shp: polygons of the smallest spatial extent at each stream site (10% of the catchment area, closest to the sampling location)  
    * forest_disturbance.shp: polygons of forest disturbance from logging, insect outbreaks, fores fire, and a general "cleared" category within the site catchments  
    * paved_roads.shp: lines of all paved roads within the site catchments
    * trails.shp: lines of all trails within the site catchments  
    * unpaved_roads.shp: lines of all unpaved roads (including ATV trails) within the site catchments  
* [surber_sampling.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/surber_sampling.csv): Number of surber samples collected from each stream site.
* [tn.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/tn.csv): Dissolved nitrogen data from filtered water samples collected at each stream site, measured with a DOC/TDN analyzer.
* [water_chemistry.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/water_chemistry.csv): Measurements of ph, water temperature, electrical conductivity, total dissolved solids, alkalinity, and turbidity from each stream site.

#### Not in repository:  
* invert_coeficients.csv: Coefficients for converting benthic invertebrate length to mass using the power law allometric equation (Burgherr and Meyer, 1997). We used coefficients from (Benke et al., 1999) for these equations, using the “all insect” category for orders where no other coefficients were available (i.e., *collembola*, *oligochaeta*, *gastropoda*, *hirudinea*, *acarina*, *neuropteran*, *lepidoptera*, and *bivalvia*). Coefficients can be found through the link in the "Data sources" section above.
* invert_functional_groups.csv: Functional groups assigned to each taxa using data from the USGS benthic invertebrate database (see "Data sources: section above)
* invert_traits_usgs.csv: File containing mean length values for each benthic invertebrate taxa from the USGS benthic invertebrate database (see "Data sources: section above). These data were used in the power law allometric equation (Burgherr and Meyer, 1997).
* tss_filters.csv: Mass of total suspended solids measured from water samples collected at each stream site. Note that these measurements were below the instrument detection limit and were not included in our analyses.
* periphyton_afdm.csv: Ash free dry mass (AFDM) of periphyton samples collected at each stream site

#### References:  

Benke, A. C., Huryn, A. D., Smock, L. A., & Wallace, J. B. (1999). Length-Mass Relationships for Freshwater Macroinvertebrates in North America with Particular Reference to the Southeastern United States. In Source: Journal of the North American Benthological Society (Vol. 18, Issue 3).  

Burgherr, P., & Meyer, E. I. (1997). Regression analysis of linear body dimensions vs. dry mass in stream macroinvertebrates. Archiv Für Hydrobiologie, 139(1), 101–112. https://doi.org/10.1127/archiv-hydrobiol/139/1997/101  

Hauer, F. R., & Lamberti, G. A. (2007). Methods in stream ecology. Elsevier Inc.  

Ministry of Environment. (2009). The Canadian Aquatic Biomonitoring Network Field Manual. http://www.unb.ca/cri/cabin_criweb.html  

### Folder 2: scripts
Scripts used for processing/analyzing our *in situ* and geospatial data and for simulating disturbances in our meta-ecosystem model.  
*  [calculating_disturbance_metrics.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/calculating_disturbance_metrics.R): Processing geospatial data (i.e., forest distubance, roads and trails, landcover, and human impact index) to calculate metrics for each site at the catchment, riparian, and local spatial extent. Generates the following files: disturbance_data_large.csv, disturbance_data_med.csv, disturbance_data_small.csv    
*  [correlogram.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/correlogram.R): script to generate the correlogram that was used to select non-correlated variables for our general linear models.
*  [global_sensitivity_analysis.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/global_sensitivity_analysis.R): Script for running the global sensitivity analysis to determine the most important parameters in our meta-ecosystem model.  
*  [meta_ecosystem_analytical_equilibria.nb](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/meta_ecosystem_analytical_equilibria.nb): Code for generating the analytical equilibria and jacobian matrix for our meta-ecosystem model, written in Wolfram Language.
*  [simulate_disturbance.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/simulate_disturbance.R): Script for simulating disturbance in our meta-ecosystem model, generates param_simulations_stable_10000.csv  
*  [stats_models.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/stats_models.R): Script to evaluate top general linear models for each key stream response variable (i.e., benthic invertebrate biomass, EPT index, periphyton biomass, percent shredders, dissolved nitrogen, electrical conductivity, and embeddedness) to determine relationships between stream quality and terrestrial disturbance. Generates empirical_glm_results.csv  
*  [stream_data_cleaning.R](https://github.com/hfadams/meta_ecosystem_model/blob/main/scripts/stream_data_cleaning.R): Script to clean and process *in situ* empirical data collected at each stream site, generates empirical_stream_data.csv and empirical_stream_data_standard_deviations.csv  

### Folder 3: output
Key files generated by the various scripts, including statistsical analysis of the empirical data and simulations generated by the meta-ecosystem model  
* [disturbance_data_catchments.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/disturbance_data_catchment.csv): Disturbance metrics calculated from geospatial data at the catchment extent
* [disturbance_data_riparian.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/disturbance_data_riparian.csv): Disturbance metrics calculated from geospatial data at the riparian extent
* [disturbance_data_local.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/disturbance_data_local.csv): Disturbance metrics calculated from geospatial data at the local extent
* [empirical_glm_results.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/empirical_glm_results.csv): 
* [empirical_stream_data.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/empirical_stream_data.csv): Mean values from the triplicate samples collected *in situ* at each stream site  
* [empirical_stream_data_standard_deviations.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/empirical_stream_data_standard_deviations.csv): Standard deviation of the triplicate samples collected *in situ* at each stream site to be paired with "empirical_stream_data.csv") 
* [param_simulations_stable.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/output/param_simulations_stable.csv): Equilibrium states for the metaecosystem using 10000 randomly generated sets of parameters, filtered for stability. The first 10000 equilibria are selected from this file and used as the initial state for the disturbance simulations.

## Sharing and accessing the data
This project is licensed under the MIT license, please see the [MIT license web page](https://opensource.org/license/mit/) for details.

## Funding
This work was funded by an NSERC Discovery grant. We would like to thank all the institutions and authors who made their data open source and free to support our work.

## Recommended citation

## Authors
### Scripts
**Hannah Adams** - *Author* - [LinkedIn](https://www.linkedin.com/in/hannah-adams-624122219/), [GitHub](https://github.com/hfadams), [ORCiD](https://orcid.org/0000-0003-2647-8021)

### Manuscript
**Hannah Adams** - *Author* - [LinkedIn](https://www.linkedin.com/in/hannah-adams-624122219/), [GitHub](https://github.com/hfadams), [ORCiD](https://orcid.org/0000-0003-2647-8021)  
**Shawn J. Leroux** - *Co-author* - [GitHub](https://github.com/sjleroux), [ORCiD](https://orcid.org/0000-0001-9580-0294), [website](https://shawnleroux.wixsite.com/lerouxlab)

## References

