rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum","class","order","family","genus")

for ( t in taxa )
{
    inFileName <- paste( "PERMANOVApValuesTaxaVsMetadata_", t, ".txt", sep="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    myColClasses <- c("character", "numeric", "numeric")
    myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

    pdf( paste(t , "_PERMANOVAplots.pdf",sep=""))
    hist(myT[,2], breaks=20, xlim=c(0, 1.0), xlab="p-value", main=t, las=1)
    dev.off()
}
