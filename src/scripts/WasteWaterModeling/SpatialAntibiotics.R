rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")
library("rgl")

## Testing with another comment
## Load Support Packages
## Load Data
## Visualize Residuals
## Hypothesize Models
## Will have a spatial and time component.
## conc ~ space*time
## Run Simple Models

setwd("/Users/mbrown67/Google\ Drive/Urban\ Environmental\ Genomics\ Project/")

inFileName <- "MetaMasterAllinOne.csv"

myT <- read.csv(inFileName, header=TRUE, na.strings = c("n/a"))#, "ND", "<LOQ"))

## Restrict it to timepoint 1
myT <- myT[myT$Time == 1,]

## Restrict it to lane6
myT <- myT[myT$Lane == 6,]

myT <- myT[myT$Rep == 1,]

## Restrict to Sugar Creek
myT <- myT[which(grepl("MALLARD", x = myT$Sample.ID)),]

## ID and spatial info

id <- myT$Sample.ID
loc <- myT$Sample.Location
lat <- myT$Lattitude
long <- myT$Longitude

time <- myT$Timepoint
rep <- myT$Rep

seqID <- myT$SampleID
seqIndex <- myT$Index
lane <- myT$Lane
## Antibiotic Info
## Find out meanings for ND and LOQ
## This will involve imputation
## This will make rarification interesting
## How will this be involved: "A polynomial curve was used for ciprofloxacin calculation"
## antibioticMatrix <- myT[,16:25]
erta <- myT$Ertapenem..ng.L.
erta <- gsub("ND|<LOQ", "0", erta)
erta <- as.numeric(erta)
amox <- myT$Amoxicillin..ng.L.
amox <- gsub("ND|<LOQ", "0", amox)
amox <- as.numeric(amox)
cipro <- myT$Ciprofloxacin..ng.L.
cipro <- gsub("ND|<LOQ", "0", cipro)
cipro <- as.numeric(cipro)
doxy <- myT$Doxycycline..ng.L.
doxy <- gsub("ND|<LOQ", "0", doxy)
doxy <- as.numeric(doxy)
azith <- myT$Azithromycin..ng.L.
azith <- gsub("ND|<LOQ", "0", azith)
azith <- as.numeric(azith)
## clinda is missing data right now; Ignore it
clinda <- myT$Clindamycin..ng.L.
clinda <- gsub("ND|<LOQ", "0", clinda)
clinda <- as.numeric(clinda)
sulfa <- myT$Sulfamethoxazole..ng.L.
sulfa <- gsub("ND|<LOQ", "0", sulfa)
sulfa <- as.numeric(sulfa)
cepha <-  myT$Cephalexin..ng.L.
cepha <- gsub("ND|<LOQ", "0", cepha)
cepha <- as.numeric(cepha)
trimetho <- myT$Trimethoprim..ng.L.
trimetho <- gsub("ND|<LOQ", "0", trimetho)
trimetho <- as.numeric(trimetho)
levo <- myT$Levofloxacin..ng.L.
levo <- gsub("ND|<LOQ", "0", levo)
levo <- as.numeric(levo)

ABlist <- list(erta, amox, cipro, doxy, azith, clinda, sulfa, cepha, trimetho, levo)

modelList <- list()
iter <- 1
for(x in ABlist){
    ##print(x)
    modelList[[iter]] <- anova(lm(as.numeric(x) ~  lat*long))
    iter <- iter + 1
}

## plot(erta, amox)

## plot3d(lat, long, erta)
## I'm not sure if I can use the standard procedure to extending to multiple antibiotics
