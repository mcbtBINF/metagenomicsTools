rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

## Testing with another comment
## Load Support Packages
## Load Data
## Visualize Residuals
## Hypothesize Models
## Will have a spatial and time component.
## conc ~ space*time
## Run Simple Models
### taxa ~ antibioticConc + (1|timePoint)

## Model 1
## taxa ~ SummedAntibiotic Concentration

## Model 2
## taxa ~ Environmental Features

## Model 3
## diversity ~ SummedAntibiotic Concentration

## Model 4
## MDS Axis ~ SummedAntibiotic Concentration

