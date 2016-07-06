
rm(list=ls())

# setwd("C:\\Caroll_Nov_2015\\spreadsheets")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum","class","order","family","genus")

pdf("allHistograms.pdf")
par(mfrow=c(3,2))


for ( t in taxa )
{
	inFileName <- paste( "samp_cohort_pValuesTaxaVsMetadata_", t, ".txt", sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",2),"numeric", "character", rep("numeric", numCols-4))
        myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

       	### pdf( paste(t , "_KendallpValueHistogram.pdf",sep=""))

        hist(myT[,3], xlim = c(0, 1.0), main=t, xlab = "p-value")
##        dev.off()
}
dev.off()
