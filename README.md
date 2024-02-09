# Global Haemosporidian Diversity
This repository contains the scripts used in the analysis presented in the manuscript titled "Energy Input, Habitat Heterogeneity, and Host Specificity Drive Avian Haemosporidian Diversity at Continental Scales"
The four scripts included are:
- `dependent variables` calulate pylogenetic diversity metrics of parasite assemblages.
- `predictors` calculates predictors such as host richness and degree of generalism and merges it with environmental data obtained from different products.
- `varimpanalysis` performs a variable importance analysis.
- `geogammodels` performs the model building process.

Parasite data was downloaded from MalAvi (Host and Sites table.
Parasite phylogenetic hypotheses were obtained using MrBayes at CIPRES geteway, followwing the workflow as impemented by the R package "treespace", 100 mcc trees were produced by identifying 100 tree clusters from a 10000 trees sample of the mcmc runs. 

Environmental datum were obtained from several sources, including:
- Bioclimatic variables of worldclim (annual mean temperature-Bio1, temperature annual range-Bio 7, precipitation seasonality-Bio 15, and annual precipitation-Bio 12).
- Net radiation from CERES product of Nasa Earth Observations (https://neo.gsfc.nasa.gov/).
- Land cover types from MODIS (MOD12Q1).  
- Global human population density (GPW) from CIESIN.
- Actual evapotranspiration (AET) from USGS.
- Human footprint from venter et al. 2017.
- Ecosystem heterogeneity calculated with Rao´s Q, diversity metric of evi measures within a cell using RASTERDIV package workflow.
