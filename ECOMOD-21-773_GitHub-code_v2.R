## *****************************************************************************
## Comparing distribution of harbour porpoise using Generalized Additive Models 
## and hierarchical Bayesian models with Integrated Nested Laplace Approximation 
## *****************************************************************************

## Setting-up
## *****************************************************************************

## Load required packages into R
library(sp)
library(Matrix)
library(INLA)
library(inlabru)
library(mgcv)

## Load data
load("gitporpoise_v2.RData")

## *****************************************************************************

## Gridded Model
## *****************************************************************************

## Specify SPDE with default priors
matern <- inla.spde2.matern(mesh)

## Specify model components
cmp_grid <- Size ~ Intercept(1) +
  field(coordinates, model=matern) + ## Spatial field using SPDE
  DepthCentered(cov_pxl, model="linear") + ## Covariate
  SP_HU3Centered(cov_pxl, model="linear") + ## Covariate
  ChlCentered(cov_pxl, model="linear") ## Covariate

## Fit model
fit_grid <- bru(components = cmp_grid, ## Model components specified above
                data = griddata, ## Analyse gridded data
                family = "nbinomial", ## Negative Binomial distribution
                options = list(E=griddata$Effort, ## Correct for sampling effort
                               control.compute=list(dic=T, ## Calculate DIC
                                                    waic=T, ## Calculate WAIC
                                                    cpo=T, ## Calculate CPO
                                                    config=T)), ## Store internal GMRF approximations 
                domain = list(coordinates=mesh)) ## Spatial domain

## Make predictions onto pred_pxl using fitted model
pred_grid <- predict(fit_grid, pred_pxl, ~ exp(Intercept +
                                                 field +
                                                 DepthCentered +
                                                 SP_HU3Centered +
                                                 ChlCentered),
                     n.samples=1000)

## *****************************************************************************

## LGCP Model
## *****************************************************************************

## Specify model components
cmp_lgcp <- coordinates ~ Intercept(1) +
  field(coordinates, model=matern) + ## Spatial field using SPDE
  SlopeCentered(cov_pxl, model="linear") + ## Covariate
  PsndgrvsndCentered(cov_pxl, model="linear") + ## Covariate
  SP_HU3Centered(cov_pxl, model="linear") + ## Covariate
  ChlCentered(cov_pxl, model="linear") + ## Covariate
  V_ShearCentered(cov_pxl, model="linear") ## Covariate

## Fit model
fit_lgcp <- lgcp(components = cmp_lgcp, ## Model components specified above
                 data = pointdata, ## Analyse point data
                 samplers = samplers, ## Sampling took place in these locations
                 domain = list(coordinates=mesh), ## Spatial domain
                 options = list(control.compute=list(dic=T, ## Calculate DIC
                                                     waic=T, ## Calculate WAIC
                                                     cpo=T, ## Calculate CPO
                                                     config=T))) ## Store internal GMRF approximations

## Make predictions onto pred_pxl using fitted model
pred_lgcp <- predict(fit_lgcp, pred_pxl, ~ exp(Intercept +
                                                 field +
                                                 SlopeCentered +
                                                 PsndgrvsndCentered +
                                                 SP_HU3Centered +
                                                 ChlCentered +
                                                 V_ShearCentered),
                     n.samples=1000)

## *****************************************************************************

## GAM 
## *****************************************************************************

## Fit GAM
fit_GAM <-gam(Size~
                DepthCentered + 
                SP_HU3Centered +
                ChlCentered +
                SlopeCentered +
                BTCentered +
                SPEEDCentered +
                Fr_SideCentered +
                V_ShearCentered +
                s(Longitude,Latitude, bs="so", xt=list(bnd=soap_boundary)), 
              knots=soap_knots, 
              offset=offset, 
              data=gamdata, 
              family=nb(), 
              method="ML")

summary(fit_GAM)
AIC(fit_GAM)
