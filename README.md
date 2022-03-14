# gitporpoise
This supplementary material provides the code used to fit the GAM, INLA gridded, and INLA LGCP analyses discussed in 'Comparing distribution of harbour porpoise using Generalized Additive Models and hierarchical Bayesian models with Integrated Nested Laplace Approximation' submitted to Ecological Modelling.

It also contains the data and model objects used in these analyses. This is available in the file gitporpoise.Rdata

## Data

* gitporpoise$pointdata
  * SpatialPointsDataFrame containing locations of harbour porpoise observations from transect surveys. Analysed using LGCP model.
  * Includes data on:
    * coordinates - Coordinate locations of observations in a Lambert azimuthal equal-area (LAEA) projection (units=km).
    * OBJECTID - ?
    * Location - ?
    * Survey_Dat - Date of survey in DD/MM/YYYY.
    * Transect - Transect ID number.
    * Size - ?
    * Year - Year of survey in YYYY.
* gitporpoise$griddata
  * SpatialPointsDataFrame containing gridded information for harbour porpoise observations and environmental data. Analysed using GAM and gridded model.
  * Includes data on:
    * Depth
    * Slope
    * Psndgrvsnd
    * X
    * Y
    * SP_HU3
    * Strip_Widt
    * Transect
    * BT
    * Chl
    * netPP
    * PEA
    * SPEED
    * SST
    * V_Shear
    * Vert_Vel
    * Wsurf
    * Size
    * Effort_Tot
    * chl_NEODAAS
    * Fr_Dist
    * Fr_Side
    * sst_NEODAAS
    * eff
    * DepthCentered - Depth (m) (centered by subtracting mean)
    * SP_HU3Centered - Mixing (H/u^3) (centered by subtracting mean)
    * BTCentered - Bottom temperature (⁰C) (centered by subtracting mean)
    * ChlCentered - Chlorophyll-aI (mgC/m^3) (centered by subtracting mean)
    * SPEEDCentered - Current speed (m/s) (centered by subtracting mean)
    * V_ShearCentered - Vertical shear (m/s) (centered by subtracting mean)
    * Fr_SideCentered - Thermal front side (centered by subtracting mean)
    * SlopeCentered - Slope (⁰) (centered by subtracting mean)
    * netPPCentered
    * PsndgrvsndCentered - Proportion of sediment that was sand or gravelly sand (%) (centered by subtracting mean)
    * PEACentered
    * SSTCentered
    * Vert_VelCentered
    * Fr_DistCentered
    * Longitude
    * Latitude 
* gitporpoise$samplers
  * SpatialLinesDataFrame containing information on line transects, used in LGCP model.
* gitporpoise$mesh
  * INLA mesh triangulation used in LGCP and gridded models. 
* SpatialPixelsDataFrames containing covariate information at 1km grid resolution across spatial area. Including:
  * gitporpoise$pxl_DepthCentered
  * gitporpoise$pxl_SP_HU3Centered
  * gitporpoise$pxl_BTCentered
  * gitporpoise$pxl_ChlCentered
  * gitporpoise$pxl_SPEEDCentered
  * gitporpoise$pxl_V_ShearCentered
  * gitporpoise$pxl_Fr_SideCentered
  * gitporpoise$pxl_SlopeCentered
  * gitporpoise$pxl_netPPCentered
  * gitporpoise$pxl_PsndgrvsndCentered
  * gitporpoise$pxl_PEACentered
  * gitporpoise$pxl_SSTCentered
  * gitporpoise$pxl_Vert_VelCentered
  * gitporpoise$pxl_Fr_DistCentered
* gitporpoise$pred_pxl
  * SpatialPixelsDataFrame containing covariate information. Used for making predictions from LGCP and gridded models. 
* gitporpoise$soap_boundary
  * Boundary for GAM soap film smoother
* gitporpoise$soap_knots
  * Knots for GAM soap film smoother
