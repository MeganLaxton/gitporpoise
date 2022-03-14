# gitporpoise
This supplementary material provides the code used to fit the GAM, INLA gridded, and INLA LGCP analyses discussed in 'Comparing distribution of harbour porpoise using Generalized Additive Models and hierarchical Bayesian models with Integrated Nested Laplace Approximation' submitted to Ecological Modelling.

It also contains the data and model objects used in these analyses. This is available in the file gitporpoise.Rdata

## Data

* pointdata
  * SpatialPointsDataFrame containing locations of harbour porpoise observations from transect surveys. Analysed using LGCP model. Includes data on:
    * coordinates - Coordinate locations of observations in a Lambert azimuthal equal-area (LAEA) projection (units=km).
* griddata
  * SpatialPointsDataFrame containing gridded information for harbour porpoise observations. Analysed using gridded model. Includes data on:
    * coordinates - Coordinate locations of observations in a Lambert azimuthal equal-area (LAEA) projection (units=km).
    * Size - number of sightings in each grid cell. 
    * Effort - length of survey effort in each grid cell.
* gamdata
  * data.frame containing gridded information for harbour porpoise observations and environmental data. Analysed using GAM. Includes data on:     
    * Longitude - Longitude coordinate location of observations.
    * Latitude - Latitude coordinate location of observations.
    * Size - number of sightings in each grid cell. 
    * offset
    * DepthCentered - Depth (m) (centered by subtracting mean)
    * SP_HU3Centered - Mixing (H/u^3) (centered by subtracting mean)
    * ChlCentered - Chlorophyll-aI (mgC/m^3) (centered by subtracting mean)
    * SlopeCentered - Slope (⁰) (centered by subtracting mean)
    * BTCentered - Bottom temperature (⁰C) (centered by subtracting mean)
    * SPEEDCentered - Current speed (m/s) (centered by subtracting mean)
    * Fr_SideCentered - Thermal front side (centered by subtracting mean)
    * V_ShearCentered - Vertical shear (m/s) (centered by subtracting mean)
* samplers
  * SpatialLinesDataFrame containing information on line transects, used in LGCP model.
* mesh
  * INLA mesh triangulation used in LGCP and gridded models. 
* cov_pxl
  * SpatialPixelsDataFrame containing covariate information at 1km grid resolution across spatial area covered by mesh. Used in LGCP and gridded models. Includes data on:
    * SlopeCentered - Slope (⁰) (centered by subtracting mean)
    * PsndgrvsndCentered - Proportion of sediment that was sand or gravelly sand (%) (centered by subtracting mean)
    * SP_HU3Centered - Mixing (H/u^3) (centered by subtracting mean)
    * ChlCentered - Chlorophyll-aI (mgC/m^3) (centered by subtracting mean)
    * V_ShearCentered - Vertical shear (m/s) (centered by subtracting mean)
    * DepthCentered - Depth (m) (centered by subtracting mean)
* pred_pxl
  * SpatialPixelsDataFrame containing covariate information at 1km grid resolution across spatial area of interest (within inner mesh boundary). Used for making predictions from LGCP and gridded models. Includes data on:
    * SlopeCentered - Slope (⁰) (centered by subtracting mean)
    * PsndgrvsndCentered - Proportion of sediment that was sand or gravelly sand (%) (centered by subtracting mean)
    * SP_HU3Centered - Mixing (H/u^3) (centered by subtracting mean)
    * ChlCentered - Chlorophyll-aI (mgC/m^3) (centered by subtracting mean)
    * V_ShearCentered - Vertical shear (m/s) (centered by subtracting mean)
    * DepthCentered - Depth (m) (centered by subtracting mean)
* soap_boundary
  * Boundary for GAM soap film smoother
* soap_knots
  * Knots for GAM soap film smoother
