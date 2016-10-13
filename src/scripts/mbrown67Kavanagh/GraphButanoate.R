rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

## setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/rdpClassifications")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging")

filePrefix <- paste0(c("PICRUSt_", tissueKept, "_Butanoate_Metabolism_"), collapse="_")

taxa <- 3

mlm<- TRUE

divider <- 4

    inFileName <- paste("PICRUSt_", taxa, "_LogNormwithMetadata.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")

    numCols <- ncol(myT)
    numMetadataCols <- 6

    pdf( paste("PICRUSt_", taxa, filePrefix, "boxplots.pdf", sep=""))
