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
    outputs: empirical_stream_data.csv, empirical_stream_data_standard_deviations.csv
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
    requires: rescaled_data_local.csv, rescaled_data_riparian.csv, rescaled_data_catchment.csv  
    outputs: model_summary_local_averaged.svg, model_summary_riparian_averaged.svg, model_summary_catchment_averaged.svg, model_averaging_formatted.csv  

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
### General
- dplyr
- reshape2
- svglite
- stringr
- data.table
- tibble

### Plotting
- ggplot2
- gridExtra
- lattice
  
### Spatial analysis
- sf
- raster
- elevatr
- lwgeom

### Statistical analysis
- AICcmodavg
- nlme
- wiqid

# Meta-ecosystem model
- lhs
- deSolve
- rootSolve
- randomForest

## References
  - Auguie B (2017). _gridExtra: Miscellaneous Functions for "Grid" Graphics_. R package version 2.3, <https://CRAN.R-project.org/package=gridExtra>.
  - Carnell R (2022). _lhs: Latin Hypercube Samples_. R package version 1.1.6, <https://CRAN.R-project.org/package=lhs>.
  - Dowle M, Srinivasan A (2023). _data.table: Extension of `data.frame`_. R package version 1.14.8, <https://CRAN.R-project.org/package=data.table>.
  - Hijmans R (2023). _raster: Geographic Data Analysis and Modeling_. R package version 3.6-26, <https://CRAN.R-project.org/package=raster>.
  - Hollister J, Shah T, Nowosad J, Robitaille A, Beck M, Johnson M (2023). _elevatr: Access Elevation Data from Various APIs_. doi:10.5281/zenodo.8335450 <https://doi.org/10.5281/zenodo.8335450>, R package version 0.99.0, <https://github.com/jhollist/elevatr/>.
  - Karline Soetaert, Peter M.J. Herman (2009). _A Practical Guide to Ecological Modelling. Using R as a Simulation Platform _. Springer. ISBN 978-1-4020-8623-6. Karline Soetaert (2009). _rootSolve: Nonlinear root finding, equilibrium and steady-state analysis of ordinary differential equations _. R package 1.6.
  - Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010). "Solving Differential Equations in R: Package deSolve." _Journal of Statistical Software_, *33*(9), 1-25. doi:10.18637/jss.v033.i09 <https://doi.org/10.18637/jss.v033.i09>.
  - Liaw A, Wiener M (2002). "Classification and Regression by randomForest." _R News_, *2*(3), 18-22. <https://CRAN.R-project.org/doc/Rnews/>.
  - Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). "Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption." _CRAN_. <https://easystats.github.io/report/>.
  - Mazerolle MJ (2023). _AICcmodavg: Model selection and multimodel inference based on (Q)AIC(c)_. R package version 2.3.3, <https://cran.r-project.org/package=AICcmodavg>.
  - Meredith M (2022). _mcmcOutput: Functions to Store, Manipulate and Display Markov Chain Monte Carlo (MCMC) Output_. R package version 0.1.3, <https://CRAN.R-project.org/package=mcmcOutput>.
  - Meredith M (2022). _wiqid: Quick and Dirty Estimates for Wildlife Populations_. R package version 0.3.3, <https://CRAN.R-project.org/package=wiqid>.
  - Meredith M, Kruschke J (2022). _HDInterval: Highest (Posterior) Density Intervals_. R package version 0.2.4, <https://CRAN.R-project.org/package=HDInterval>.
  - Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1, <https://CRAN.R-project.org/package=tibble>.
  - Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2024). _vegan: Community Ecology Package_. R package version 2.6-6.1, <https://CRAN.R-project.org/package=vegan>.
  - Pebesma E (2023). _lwgeom: Bindings to Selected 'liblwgeom' Functions for Simple Features_. R package version 0.2-13, <https://CRAN.R-project.org/package=lwgeom>.
  - Pebesma E, Bivand R (2005). "Classes and methods for spatial data in R." _R News_, *5*(2), 9-13. <https://CRAN.R-project.org/doc/Rnews/>. Bivand R, Pebesma E, Gomez-Rubio V (2013). _Applied spatial data analysis with R, Second edition_. Springer, NY. <https://asdar-book.org/>.
  - Pebesma E, Bivand R (2023). _Spatial Data Science: With applications in R_. Chapman and Hall/CRC. doi:10.1201/9780429459016 <https://doi.org/10.1201/9780429459016>, <https://r-spatial.org/book/>. Pebesma E (2018). "Simple Features for R: Standardized Support for Spatial Vector Data." _The R Journal_, *10*(1), 439-446. doi:10.32614/RJ-2018-009 <https://doi.org/10.32614/RJ-2018-009>, <https://doi.org/10.32614/RJ-2018-009>.
  - Pinheiro J, Bates D, R Core Team (2023). _nlme: Linear and Nonlinear Mixed Effects Models_. R package version 3.1-162, <https://CRAN.R-project.org/package=nlme>. Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. doi:10.1007/b98882 <https://doi.org/10.1007/b98882>.
  - R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
  - Sarkar D (2008). _Lattice: Multivariate Data Visualization with R_. Springer, New York. ISBN 978-0-387-75968-5, <http://lmdvr.r-forge.r-project.org>.
  - Simpson G (2022). _permute: Functions for Generating Restricted Permutations of Data_. R package version 0.9-7, <https://CRAN.R-project.org/package=permute>.
  - Wickham H (2007). "Reshaping Data with the reshape Package." _Journal of Statistical Software_, *21*(12), 1-20. <http://www.jstatsoft.org/v21/i12/>.
  - Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
  - Wickham H (2022). _stringr: Simple, Consistent Wrappers for Common String Operations_. R package version 1.5.0, <https://CRAN.R-project.org/package=stringr>.
  - Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.3, <https://CRAN.R-project.org/package=dplyr>.
  - Wickham H, Henry L, Pedersen T, Luciani T, Decorde M, Lise V (2023). _svglite: An 'SVG' Graphics Device_. R pac> ge version 2.1.2, <https://CRAN.R-project.org/package=svglite>.
