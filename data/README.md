# Empirical data files
*In situ* data collected from stream sites and geospatial data generated from satellite imagery in QGIS, used to find relationships between stream characteristics and terrestrial disturbances across spatial extents.  

## Folder directory  
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
* [tss_filters.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/tss_filters.csv): Mass of total suspended solids measured from water samples collected at each stream site. Note that these measurements were below the instrument detection limit and were not included in our analyses.
* [water_chemistry.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/water_chemistry.csv): Measurements of ph, water temperature, electrical conductivity, total dissolved solids, alkalinity, and turbidity from each stream site.
* [periphyton_afdm.csv](https://github.com/hfadams/meta_ecosystem_model/blob/main/data/periphyton_afdm.csv): Ash free dry mass (AFDM) of periphyton samples collected at each stream site

#### Not in repository:  
* invert_coeficients.csv: Coefficients for converting benthic invertebrate length to mass using the power law allometric equation (Burgherr and Meyer, 1997). We used coefficients from (Benke et al., 1999) for these equations, using the “all insect” category for orders where no other coefficients were available (i.e., *collembola*, *oligochaeta*, *gastropoda*, *hirudinea*, *acarina*, *neuropteran*, *lepidoptera*, and *bivalvia*). Coefficients can be found through the link in the "Data sources" section above.
* invert_functional_groups.csv: Functional groups assigned to each taxa using data from the USGS benthic invertebrate database (see "Data sources: section above)
* invert_traits_usgs.csv: File containing mean length values for each benthic invertebrate taxa from the USGS benthic invertebrate database (see "Data sources: section above). These data were used in the power law allometric equation (Burgherr and Meyer, 1997).

#### References:  

Benke, A. C., Huryn, A. D., Smock, L. A., & Wallace, J. B. (1999). Length-Mass Relationships for Freshwater Macroinvertebrates in North America with Particular Reference to the Southeastern United States. In Source: Journal of the North American Benthological Society (Vol. 18, Issue 3).  

Burgherr, P., & Meyer, E. I. (1997). Regression analysis of linear body dimensions vs. dry mass in stream macroinvertebrates. Archiv Für Hydrobiologie, 139(1), 101–112. https://doi.org/10.1127/archiv-hydrobiol/139/1997/101  

Hauer, F. R., & Lamberti, G. A. (2007). Methods in stream ecology. Elsevier Inc.  

Ministry of Environment. (2009). The Canadian Aquatic Biomonitoring Network Field Manual. http://www.unb.ca/cri/cabin_criweb.html  

## Metadata for variables in each file:  

### benthic_invertebrates.csv
|       Variable       |       description           |
