rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum", "class", "order", "family", "genus")

for ( t in taxa )
{
 	inFileName <- paste(t, "_SparseThreeCol.txt", sep = "")
	myT <- read.table(inFileName, header=FALSE, sep = "\t")
	myColClasses <- c(rep("character", 2), "numeric")
	myT <- read.table(inFileName, header = FALSE, sep = "\t", colClasses = myColClasses)
        colnames(myT) <- c("id", "taxa", "count")
        myT$id <- unlist(lapply(strsplit(myT$id, split="\\."), "[[", 1))
        myT$read <- unlist(lapply(strsplit(myT$id, split="_"), "[[", 1))
	myT <- myT[ myT$read == "r1", ]

	names <- vector()
	namesA <- vector()
	namesB <- vector()
	pValueLinear <- vector()
	kendallP <- vector()
	rVals <- vector()
        pValueTestvar <- vector()
        pValueCohort <- vector()
        iccCohort <- vector()

	index <- 1

        #Reduce the number of columns by one to account for the cohort column
	for( i in c(3, taxCol : ( ncol(myT) - 1 ) ) )
	{
		if( sum( myT[,i] > 0 ) > nrow(myT) /4 )
                    {

			 for ( j in 5:(taxCol - 1))
			 {
			 	namesA[index] <- names(myT)[i]
			 	namesB[index] <- names(myT)[j]
			 }
		}
	}

write.table( file = paste( "HighAbundance_", t, ".txt", sep = ""), dFrame, row.names = FALSE, sep = "\t")
}

