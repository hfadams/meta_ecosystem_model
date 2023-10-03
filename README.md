# Meta-ecosystem model

The scripts in this repository are used to 1) clean and analyze empirical data collected from streams in Gros Morne National Park and Terra Nova National Park in Newfoundland, Canada and 2) simulate terrestrial disturbance in a meta-ecosystem model of these stream-riparian ecosystems. These files are inteded to make my analyses as transparent and reproducible as possible, and have the stream data openly accessible for further analysis. A thorough description of the collection methods can be found in the appendix of the associated thesis found in the MUN archives: [INSERT LINK WHEN AVAILABLE]. 

The land class data are not publically available, but Forest Resource Inventory shapefiles can be requested from the Newfoundland government, or can be digitized from satellite imagery.

## Data sources  

1) Field collection July-August 2022  
Water quality and stream measurements, found in "data" folder  

2) United States Geological Survey (USGS) invertebrate database  
Vieira, N.K.M. et al. (2016) ‘A Database of Lotic Invertebrate Traits for North America’, U.S. Geological Survey Data Series, 187, pp. 1–15.  
[USGS Database](https://doi.org/10.3133/ds187)

4) CanElevation 5m resolution digital elevation model  
Government of Canada (2022) High Resolution Digital Elevation Model (HRDEM) - CanElevation Series.
[CanElevation HRDEM](https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995) 

6) Data agreement with Gros Morne National Park and Terra Nova national park  
Gros Morne contact information:  
Terra Nova contact information:  
    
7) Data agreement with the provincial government  
FRI data request: 

## Methods  
### In situ data  
I collected the in situ data in triplicate at each stream site and took the mean of the three samples. Refer to the main thesis text and appendix A for detailed methods.

### Model simulations
After developing a riparian-stream meta-ecosystem model, I used Mathematica to solve for all possible analytical equilibria of the model. I then selected the equilibrium that was locally stable and feasible, and parameterized the model by generating 10,000 random parameter combinations, each within a range from 0-10 (or 0-1 if a proportion). From these I selected the first 1,000 equilibria taht were feasible, locally stable, and where the benthic invertebrate biomass was greater than periphyton biomass. I used these simulations as the "undisturbed" meta-ecosystem to which I created "terrestrial disturbances" by increasing key parameters to simulate tree removal and increased erosion. Refer to the main thesis text and Appendix B for further details.

### Quality assurance  
I removed all data below the method detection limit of each *in situ* measurement and lab analysis before statistical analysis. I ensured that there was no correlation between predictor variables in the empirical dataset before developing the general linear models. 

I performed a global sensitivity analysis on each trophic level and productivity metric in the meta-ecosystem model to identify parameters creating the most uncertainty in the model.  

## software and packages
All data processing and analyses for this project were implemented using R (ver. 4.2.2), Mathematica (ver. 13.2.1), and QGIS (ver. 3.26.3).

## repository directory
### Folder 1: original data

### Folder 2: scripts

### Folder 3: output

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
