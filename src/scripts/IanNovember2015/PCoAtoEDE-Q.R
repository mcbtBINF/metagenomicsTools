rm(list=ls())

# setwd("C:\\Caroll_Nov_2015\\spreadsheets")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- "phylum"#,"class","order","family","genus")

unifrac.pc <- read.table("NEWpcs_1-26-2016.txt", sep="\t", header=FALSE)

selectedSamplesPC <- c(65, 94, 43, 58, 84, 31, 10, 15, 4, 26, 68, 92, 23, 60, 83, 59, 82, 38,
                       45, 74, 48, 20, 3, 2, 13, 34, 57, 87, 100, 32, 30, 56, 70, 75, 44, 12,
                       93, 86, 33, 5, 53, 49, 71, 42, 6, 41, 16, 80, 78, 39, 24, 76, 8, 69, 47,
                       55, 19, 25, 37, 99, 67, 27, 51, 9, 14, 91, 54, 79, 81, 73, 66, 22, 85, 11,
                       1, 29, 96, 35, 98, 50, 18, 95, 90, 52, 40, 17, 72, 63, 21, 28, 64)


    	inFileName <- paste( "pivoted_", taxa, "asColumnsLogNormalPlusMetadata.txt" , sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",2),"numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
	myT <- myT[ myT$read== "r1" & ! is.na(myT$AGE) , ]
        myT <- myT[as.numeric(sub("r1_", "", myT$id)) %in% selectedSamplesPC,]

      	pdf( paste(taxa , "samp_cohort_correlationwithEATINGPlots.pdf",sep=""))

        plot(unifrac.pc[,1], myT$EATING)
        dev.off()
