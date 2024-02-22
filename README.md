# Energy input, habitat heterogeneity, and host specificity on avian haemosporidian diversity at continental scales

<https://doi.org/10.5061/dryad.hx3ffbgkg>
This repository contains the scripts necessary to reproduce the analysis presented in the homonymous paper published in Proceedings B.

## Description of the data and file structure

This dataset is composed of 5 R scripts (*.R), one script contains functions and the other four contain the workflow to reproduce the analysis, the relationships between data and scripts are described below.

*   Note: The first two scripts load the files from a folder called "data", and the resulting files are saved in a folder called "out". In addition, all scripts load functions from a subfolder called "functions" within a folder called "scripts". Make sure to modify these paths in the scripts or create these folders with the corresponding files (input) on your workspace.

### Details for  `1diversity metrics.R`

#### Description:

Calculates phylogenetic diversity metrics of parasite assemblages: taxonomic diversity, residual phylogenetic diversity, and phylogenetic structure.

#### Input files:

*   `plasmodiumPAM.csv` Presence - absence matrix for Plasmodium genus
*   `plas100trees` multyphylo object with 100 phylogenetic trees of Plasmodium lineages
*   `haemoproteusPAM.csv` Presence - absence matrix for Haemoproteus genus
*   `hae100trees` multyphylo object with 100 phylogenetic trees of Haemoproteus lineages
*   `leucocytozoonPAM.csv` Presence - absence matrix for Leucocytozoon genus
*   `leuc100trees` multyphylo object with 100 phylogenetic trees of Leucocytozoon lineages

#### Output files:

*   `plassr.csv` Taxonomic richness for Plasmodium assemblages
*   `plaspsv100.csv` 100 Measurements of phylogenetic structure of Plasmodium assemblages
*   `plasrpd100.csv` 100 Measurements of residual phylogenetic diversity of Plasmodium assemblages
*   `haesr.csv` Taxonomic richness for Haemoproteus assemblages
*   `haepsv100.csv` 100 Measurements of phylogenetic structure of Haemoproteus assemblages
*   `haerpd100.csv` 100 Measurements of residual phylogenetic diversity of Haemoproteus assemblages
*   `leusr.csv` Taxonomic richness for Leucocytozoon assemblages
*   `leupsv100.csv` 100 Measurements of phylogenetic structure of Leucocytozoon assemblages
*   `leurpd100.csv` 100 Measurements of residual phylogenetic diversity of Leucocytozoon assemblages

### Details for  `2degree of generalism.R`

#### Description:

Measures the degree of generalism of parasite assemblages and merges with a dataset with all other studied predictors.

#### Input files:

*   `MXCLADECREDBIRDTREE.tre` maximum clade credibility tree of host species

*   `plasmodiumPAM.csv` Presence - absence matrix for Plasmodium genus

*   `haemoproteusPAM.csv` Presence - absence matrix for Haemoproteus genus

*   `leucocytozoonPAM.csv` Presence - absence matrix for Leucocytozoon genus

*   `hostxlineages.csv`   matrix for interactions between parasite lineages and host species based on Malavi records

*   `malavisubset.csv` Dataframe with malavi records for the studied regions, it includes columns such as "Lineage_Name","parasiteGenus" and coordinates "long." and "lat."

*   `regionsgrid.shp` Grid (0.25x0.25) for the extent of the studied zoogeographic realms as classified by (Holt 2014).

*   `envpredictors.csv` Dataframe with other studied predictors, dataset includes "x" and "y" as centroid coordinates (Longitude and latitude) for each cell. "aet" is the mean actual evapotranspirations obtained from USGS,"Ec-Het" is the Ecosystem heterogeneity calculated using rasterdiv package. "Humanpopdens" is the mean human population density from the GPW dataset provided by CIESIN, "Human_footprint" is the mean Human footprint obtained from data published by Venter et al. 2017. "Prec.seas." is the precipitation seasonality-Bio 15 obtained from worldclim, "rad" is the Net radiation obtained from CERES product, "shannon_diversity" is the landscape heterogeneity measured by the Shannon index calculated for diversity of ladncovers obtained from MOD12Q1 product from MODIS, "trange" is the temperature annual range-Bio 7 from worldclim, "Temperature" was obtained from annual mean temperature-Bio1 from worldclim, "precipitation" is the mean annual precipitation obtained from worlclim, "Host_richness" is the host species richness calculated as the number of host species present in each cell, this was calculated using data derivated from birdlife. With the exception of spatial coordinates, all data in this dataframe was scaled.

#### Output files:

*   `plaspredictors.csv` Dataframe generated from merging degree of generalism and environmental predictors dataframe for Plasmodium assemblages

*   `haepredictors.csv` Dataframe generated from merging degree of generalism and environmental predictors dataframe for Haemoproteus assemblages

*   `leupredictors.csv` Dataframe generated from merging degree of generalism and environmental predictors dataframe for Leucocytozoon assemblages

### Details for  `3varimpanalysis.R`

#### Description:

Performs variable importance analysis for measuring the relative importance of each studied predictor for explaining variation in parasite diversity metrics.

#### Input data:

*   Diversity metrics for each genus produced in `1diversity metrics.R` (`plaspsv100.csv`,
    `plasrpd100.csv`, `plassr.csv`, `haepsv100.csv`,
    `haerpd100.csv`, `jaesr.csv`, `leupsv100.csv`,
    `leurpd100.csv`, `leusr.csv`)
*   Dataframes with predictors produced in `2degree of generalism.R` (`plaspredictors.csv`, `haepredictors.csv`, `leupredictors.csv`)

#### Output files:

*   `todelsr.cvs` Predictors to delete for the model building process for species richness of the parasite assemblages of the three haemosporidian genera
*   `todelrpd.cvs` Predictors to delete for the model building process for residual phylogentic diversity of the parasite assemblages of the three haemosporidian genera
*   `todelpsv` Predictors to delete for the model building process for phylogenetic structure of the parasite assemblages of the three haemosporidian genera

### Details for  `4geogammodels.R`

#### Description:

Performs the model building process

#### Input data:

*   Diversity metrics for each genus produced in `1diversity metrics.R` (`plaspsv100.csv`,
    `plasrpd100.csv`, `plassr.csv`, `haepsv100.csv`,
    `haerpd100.csv`, `jaesr.csv`, `leupsv100.csv`,
    `leurpd100.csv`, `leusr.csv`)
*   Dataframes with predictors produced in `2degree of generalism.R` (`plaspredictors.csv`, `haepredictors.csv`, `leupredictors.csv`)
*   Predictors to be excluded produced in `3varimpanalysis.R` (`todelsr.cvs`, `todelrpd.cvs`, `todelpsv.cvs`)

## Sharing/Access information

All data used in this analysis was derived, transformed or adapted from the following sources:

*   MALAVI (<http://130.235.244.92/Malavi/>)
*   Worldclim (<https://www.worldclim.org/>)
*   CERES product of Nasa Earth Observations (<https://neo.gsfc.nasa.gov/>)
*   MODIS product MOD12Q1 (<https://lpdaac.usgs.gov/products/mcd12q1v061/>)
*   SEDAC-CIESIN (<https://sedac.ciesin.columbia.edu/data/collection/gpw-v4>)
*   USGS (<https://data.usgs.gov/datacatalog/data/USGS:5e39c336e4b0a79317e15d7a>)
*   Venter et al. 2017. (<https://datadryad.org/stash/dataset/doi:10.5061/dryad.052q5>)
*   Birdlife (<https://www.birdlife.org/>)


## Code/Software

All code was produced using R software version 4.2.2 in RStudio version 2023.9.1.494

### libraries for each script:

#### `1diversity metrics.R`

*   picante, dplyr, doParallel, parallel, foreach, mgcv

#### `2degree of generalism.R`

*   ape, picante, dplyr, reshape2, rgdal, raster

#### `3varimpanalysis.R`

*   dplyr, parallel, mgcv, caret, doParallel

#### `4geogammodels.R`

*   doParallel, dplyr, parallel, geoGAM, mgcv, caret, ggplot2

