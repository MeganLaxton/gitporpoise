## *****************************************************************************
## Comparing distribution of harbour porpoise using Generalized Additive Models 
## and hierarchical Bayesian models with Integrated Nested Laplace Approximation 
## *****************************************************************************

## This file formats data which is then saved as an .Rdata object to be loaded 
## in to GitHub alongside model code. Required objects are:
## - Point Data
## - Gridded Data (also used in GAM)
## - Samplers transects
## - Mesh
## - Covariate SpatialPixelsDataFrame
## - Prediction SpatialPixelsDataFrame
## - Boundary for GAM soap film smoother
## - Knots for GAM soap film smoother

## Set working directory
setwd("C:/Users/2191979L/OneDrive - University of Glasgow/Research_Projects/HP_INLAvsGAMs")

## Load required packages into R
library(sp)
library(Matrix)
library(INLA)
library(rgeos)
library(raster)
library(fields)
library(mgcv)

## Load Data
load("porpoise.RData")
griddata = porpoise$griddata
mesh = porpoise$mesh
pointdata = porpoise$points
ECoast1<-read.csv("preddata_NEODAAS.csv", header=TRUE)
coast<-read.csv("coastline_for_soapfilm_smooth.csv", head=TRUE)

## Transform point data to match projection of gridded data
pointdata <- spTransform(pointdata, CRSobj = griddata@proj4string) 

## Extract mesh outer boundary
inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
          calling inla.mesh2sp()."))
  }
  
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}
sp_mesh <- inla.mesh2sp(mesh=mesh)
mesh_boundary <- gUnaryUnion(sp_mesh$triangles)
mesh_boundary@proj4string <- mesh$crs

## Crop samplers transects to be within mesh outer boundary
samplers <- crop(porpoise$samplers, mesh_boundary)

## Center covariates by subtracting mean values
griddata$SlopeCentered = griddata$Slope - mean(griddata$Slope)
griddata$PsndgrvsndCentered = griddata$Psndgrvsnd - mean(griddata$Psndgrvsnd)
griddata$SP_HU3Centered = griddata$SP_HU3 - mean(griddata$SP_HU3)
griddata$ChlCentered = griddata$Chl - mean(griddata$Chl)
griddata$V_ShearCentered = griddata$V_Shear - mean(griddata$V_Shear)
griddata$DepthCentered = griddata$Depth - mean(griddata$Depth)

## Construct a raster covering area with 1km grid cells
r <- raster(ext=extent(-282,-68,5938,6228),
            crs=griddata@proj4string,
            ncol=214, nrow=290)

## Extract cell locations from raster
xy <- data.frame(xyFromCell(r, 1:ncell(r)))

## For each covariate, rasterise and interpolate using a thin plate spline model

### *****
### Depth 
### *****
r_DepthCentered <- rasterize(griddata@coords, r, griddata$DepthCentered, fun=mean)
r_DepthCentered@crs <- griddata@proj4string
#### Remove any NA values
v_DepthCentered <- getValues(r_DepthCentered)
i_DepthCentered <- !is.na(v_DepthCentered)
#### Extract values
xy_DepthCentered <- xy[i_DepthCentered,]
v_DepthCentered <- v_DepthCentered[i_DepthCentered]
#### Thin plate spline model
tps_DepthCentered <- Tps(xy_DepthCentered, v_DepthCentered)
r_DepthCentered <- interpolate(r_DepthCentered, tps_DepthCentered)
#### Convert to SpatialPixelsDataFrame
pxl_DepthCentered <- as(r_DepthCentered, "SpatialPixelsDataFrame")
names(pxl_DepthCentered)[1] <- "DepthCentered"

### ******
### SP_HU3
### ******
r_SP_HU3Centered <- rasterize(griddata@coords, r, griddata$SP_HU3Centered, fun=mean)
r_SP_HU3Centered@crs <- griddata@proj4string
#### Remove any NA values
v_SP_HU3Centered <- getValues(r_SP_HU3Centered)
i_SP_HU3Centered <- !is.na(v_SP_HU3Centered)
#### Extract values
xy_SP_HU3Centered <- xy[i_SP_HU3Centered,]
v_SP_HU3Centered <- v_SP_HU3Centered[i_SP_HU3Centered]
#### Thin plate spline model
tps_SP_HU3Centered <- Tps(xy_SP_HU3Centered, v_SP_HU3Centered)
r_SP_HU3Centered <- interpolate(r_SP_HU3Centered, tps_SP_HU3Centered)
#### Convert to SpatialPixelsDataFrame
pxl_SP_HU3Centered <- as(r_SP_HU3Centered, "SpatialPixelsDataFrame")
names(pxl_SP_HU3Centered)[1] <- "SP_HU3Centered"

### ***
### Chl
### ***
r_ChlCentered <- rasterize(griddata@coords, r, griddata$ChlCentered, fun=mean)
r_ChlCentered@crs <- griddata@proj4string
#### Remove any NA values
v_ChlCentered <- getValues(r_ChlCentered)
i_ChlCentered <- !is.na(v_ChlCentered)
#### Extract values
xy_ChlCentered <- xy[i_ChlCentered,]
v_ChlCentered <- v_ChlCentered[i_ChlCentered]
#### Thin plate spline model
tps_ChlCentered <- Tps(xy_ChlCentered, v_ChlCentered)
r_ChlCentered <- interpolate(r_ChlCentered, tps_ChlCentered)
#### Convert to SpatialPixelsDataFrame
pxl_ChlCentered <- as(r_ChlCentered, "SpatialPixelsDataFrame")
names(pxl_ChlCentered)[1] <- "ChlCentered"

### *******
### V_Shear
### *******
r_V_ShearCentered <- rasterize(griddata@coords, r, griddata$V_ShearCentered, fun=mean)
r_V_ShearCentered@crs <- griddata@proj4string
#### Remove any NA values
v_V_ShearCentered <- getValues(r_V_ShearCentered)
i_V_ShearCentered <- !is.na(v_V_ShearCentered)
#### Extract values
xy_V_ShearCentered <- xy[i_V_ShearCentered,]
v_V_ShearCentered <- v_V_ShearCentered[i_V_ShearCentered]
#### Thin plate spline model
tps_V_ShearCentered <- Tps(xy_V_ShearCentered, v_V_ShearCentered)
r_V_ShearCentered <- interpolate(r_V_ShearCentered, tps_V_ShearCentered)
#### Convert to SpatialPixelsDataFrame
pxl_V_ShearCentered <- as(r_V_ShearCentered, "SpatialPixelsDataFrame")
names(pxl_V_ShearCentered)[1] <- "V_ShearCentered"

### *****
### Slope
### *****
r_SlopeCentered <- rasterize(griddata@coords, r, griddata$SlopeCentered, fun=mean)
r_SlopeCentered@crs <- griddata@proj4string
#### Remove any NA values
v_SlopeCentered <- getValues(r_SlopeCentered)
i_SlopeCentered <- !is.na(v_SlopeCentered)
#### Extract values
xy_SlopeCentered <- xy[i_SlopeCentered,]
v_SlopeCentered <- v_SlopeCentered[i_SlopeCentered]
#### Thin plate spline model
tps_SlopeCentered <- Tps(xy_SlopeCentered, v_SlopeCentered)
r_SlopeCentered <- interpolate(r_SlopeCentered, tps_SlopeCentered)
#### Convert to SpatialPixelsDataFrame
pxl_SlopeCentered <- as(r_SlopeCentered, "SpatialPixelsDataFrame")
names(pxl_SlopeCentered)[1] <- "SlopeCentered"

### **********
### Psndgrvsnd
### **********
r_PsndgrvsndCentered <- rasterize(griddata@coords, r, griddata$PsndgrvsndCentered, fun=mean)
r_PsndgrvsndCentered@crs <- griddata@proj4string
#### Remove any NA values
v_PsndgrvsndCentered <- getValues(r_PsndgrvsndCentered)
i_PsndgrvsndCentered <- !is.na(v_PsndgrvsndCentered)
#### Extract values
xy_PsndgrvsndCentered <- xy[i_PsndgrvsndCentered,]
v_PsndgrvsndCentered <- v_PsndgrvsndCentered[i_PsndgrvsndCentered]
#### Thin plate spline model
tps_PsndgrvsndCentered <- Tps(xy_PsndgrvsndCentered, v_PsndgrvsndCentered)
r_PsndgrvsndCentered <- interpolate(r_PsndgrvsndCentered, tps_PsndgrvsndCentered)
#### Convert to SpatialPixelsDataFrame
pxl_PsndgrvsndCentered <- as(r_PsndgrvsndCentered, "SpatialPixelsDataFrame")
names(pxl_PsndgrvsndCentered)[1] <- "PsndgrvsndCentered"

## Approximate inner mesh boundary
### Load coastline data and convert to a spatial object
coastline <- read.csv("Coastline_updated_for_INLAbru.csv", header=TRUE)
coordinates(coastline) <- c("Long","Lat")
proj4string(coastline) <- "+proj=longlat"
coastline <- spTransform(coastline, griddata@proj4string)
### Create mesh boundary segment using coastline data
bnd.coast <- inla.mesh.segment(coordinates(coastline), is.bnd=TRUE)
### Create a nonconvex hull using gridded data
bnd1 <- inla.nonconvex.hull(griddata, 9)
### Create a SpatialPolygon for inner mesh boundary
bnd_inner <- SpatialPolygons(Srl=list(
  Polygons(srl=list(
    Polygon(
      coords=matrix(c(bnd1$loc[,1],
                      bnd1$loc[,2]),
                    ncol=2,
                    byrow=FALSE,
                    dimnames=list(NULL,
                                  c("x","y"))))),
    ID="a")), proj4string=griddata@proj4string)
### Create a SpatialPolygon for coastline
coastline_p  <- SpatialPolygons(Srl=list(
  Polygons(srl=list(
    Polygon(
      coords=matrix(c(bnd.coast$loc[,1],
                      bnd.coast$loc[,2]),
                    ncol=2,
                    byrow=FALSE,
                    dimnames=list(NULL,
                                  c("x", "y"))))),
    ID="a")), proj4string=griddata@proj4string)
### Remove areas from inner boundary that coastline intersects into
bnd_inner <- gDifference(bnd_inner, coastline_p) 

## Create a SpatialPixels object on which to make predictions
## Should contain all covariate information
cov_pxl <- cbind(pxl_SlopeCentered, pxl_PsndgrvsndCentered)
cov_pxl <- cbind(cov_pxl, pxl_SP_HU3Centered)
cov_pxl <- cbind(cov_pxl, pxl_ChlCentered)
cov_pxl <- cbind(cov_pxl, pxl_V_ShearCentered)
cov_pxl <- cbind(cov_pxl, pxl_DepthCentered)

pred_pxl <- cov_pxl[bnd_inner,]

## GAM *************************************************************************
porp<-read.csv("sightings_NEODAAS.csv", header=TRUE)                                 #\\Data - no outliers

porp$XMEAN.J <- jitter(porp$X, amount=1)
porp$YMEAN.J <- jitter(porp$Y, amount=1)
#porp$Slope.sqrt <- sqrt(porp$Slope)
porp$offset2<-log10(porp$Effort*(porp$Strip_Widt/1000))   ##convert strip width into km then multiply - no need to use log

porp2<-porp
porp2[19]<-NULL

ECoast1<-read.csv("preddata_NEODAAS.csv", header=TRUE)
library(fields)

## make soap film smoother
coast<-read.csv("coastline_for_soapfilm_smooth.csv", head=TRUE)
## create the boundary for the soap film
bound<-list(list(Longitude=coast[,1], Latitude=coast[,2], f=rep(0, nrow(coast))))

# to use a normal soap film smoother use this code:
# set the internat knots for the soap film
N <- 10
gx <- seq(min(as.data.frame(bound[[1]][1])), max(as.data.frame(bound[[1]][1])), len = N)
gy <- seq(min(as.data.frame(bound[[1]][2])), max(as.data.frame(bound[[1]][2])), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("Longitude","Latitude")
knots <- gp[with(gp, inSide(bound, Longitude, Latitude)), ]
names(knots) <- c("Longitude", "Latitude")
names(bound[[1]]) <- c("Longitude", "Latitude", "f")

## N=10   some points are too close to boundary so need to adjust
knots[1,2]<-56.13
knots[14,2]<-57.62
knots[23,1]<--3.4

## GAM *************************************************************************


## rename files to save
soap_boundary <- bound 
soap_knots <- knots

## Only include necessary info in data...
pointdata <- pointdata[,0]

griddata$Effort <- griddata$eff
griddata <- griddata[,c("Size","Effort")]

gamdata <- porp2
gamdata$offset <- gamdata$offset2

## Center covariates by subtracting the mean
gamdata$DepthCentered <- gamdata$Depth - mean(gamdata$Depth)
gamdata$SP_HU3Centered <- gamdata$SP_HU3 - mean(gamdata$SP_HU3)
gamdata$ChlCentered <- gamdata$Chl - mean(gamdata$Chl)
gamdata$SlopeCentered <- gamdata$Slope - mean(gamdata$Slope)
gamdata$BTCentered <- gamdata$BT - mean(gamdata$BT)
gamdata$SPEEDCentered <- gamdata$SPEED - mean(gamdata$SPEED)
gamdata$Fr_SideCentered <- gamdata$Fr_Side - mean(gamdata$Fr_Side)
gamdata$V_ShearCentered <- gamdata$V_Shear - mean(gamdata$V_Shear)

gamdata <- gamdata[,c("Longitude","Latitude","Size","offset","DepthCentered","SP_HU3Centered","ChlCentered","SlopeCentered","BTCentered","SPEEDCentered","Fr_SideCentered","V_ShearCentered")]

save(pointdata, ## Data used in LGCP model
     griddata, ## Data used in gridded model
     gamdata, ## Data used in GAM
     samplers, ## Samplers transects
     mesh, ## Mesh
     cov_pxl, ## Covariate SpatialPixelsDataFrame
     pred_pxl, ## Prediction SpatialPixelsDataFrame
     soap_boundary, ## Boundary for GAM soap film smoother
     soap_knots, ## Knots for GAM soap film smoother
     file = "C:/Users/2191979L/OneDrive - University of Glasgow/Research_Projects/HP_INLAvsGAMs/GitHub/gitporpoise_v2.RData")