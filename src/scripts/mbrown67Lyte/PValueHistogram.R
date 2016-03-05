rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")

for (taxa in taxaLevels)
{
    inFileName <- paste("pValuesForTaxa_bug_condition_cage_", taxa, ".txt",sep="")
    myT <- read.csv(inFileName, header=TRUE, sep="")
    pdf( paste(taxa, "ConditionPValueHistograms_.pdf", sep=""))
    hist(myT$pValuesAcuteChronic, main=taxa)
    dev.off()
}
