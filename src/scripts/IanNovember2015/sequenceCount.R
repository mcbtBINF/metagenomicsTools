rm(list=ls())

# setwd("C:\\Caroll_Nov_2015\\spreadsheets")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum","class","order","family","genus")

## Proceeding forward with the 91 samples from the PC analysis.
selectedSamplesPC <- c(65, 94, 43, 58, 84, 31, 10, 15, 4, 26, 68, 92, 23, 60, 83, 59, 82, 38,
                       45, 74, 48, 20, 3, 2, 13, 34, 57, 87, 100, 32, 30, 56, 70, 75, 44, 12,
                       93, 86, 33, 5, 53, 49, 71, 42, 6, 41, 16, 80, 78, 39, 24, 76, 8, 69, 47,
                       55, 19, 25, 37, 99, 67, 27, 51, 9, 14, 91, 54, 79, 81, 73, 66, 22, 85, 11,
                       1, 29, 96, 35, 98, 50, 18, 95, 90, 52, 40, 17, 72, 63, 21, 28, 64)
iter <- 1
for ( t in taxa )
{

        inFileName <- paste( t, "_SparseThreeCol.txt", sep="")
	myT <-read.table(inFileName,header=FALSE,sep="\t")
	numCols <- ncol(myT)
        myT <- myT[grepl("^r1_", myT[,1]),]
        rowNums <-  as.numeric(sapply(strsplit(as.character(sapply(strsplit(as.character(myT[,1]), split="\\."), function(x) x[1])), "r1_"), function(y) y[2]))
        myT[,1] <- rowNums
        if (iter == 1) {
            fullAgg <- aggregate(myT[,3]~myT[,1], data = myT, FUN=sum)
        }
        myT <- myT[myT[,1] %in% selectedSamplesPC,]
        print(sum(myT[,3])/length(selectedSamplesPC))
        myTAgg <- aggregate(myT[,3]~myT[,1], data=myT, FUN=sum)
        print(range(myTAgg[,2]))
        iter <- iter + 1
                                        #        selRows <- "r1_" %in% myT[,1]

                                        #	myT <- myT["r1_",]
##        myT <- myT[as.numeric(sub("r1_", "", myT$id)) %in% selectedSamplesPC,]

}

