rm(list=ls())

# setwd("C:\\Caroll_Nov_2015\\spreadsheets")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum","class","order","family","genus")

nPercent <- 0.01

### Samples selected come from the 91 in the MDS plots
### These represent out of the 100 beginning participants, those that 1) completed the various psych surveys 2) submitted a fecal sample and 3) the reads mapped (at the phyla threshold) met a cutoff of 2500. In fact, the 3 sequences excluded for depth reasons all had less than 100 reads mapping at the phyla level.
selectedSamplesPC <- c(65, 94, 43, 58, 84, 31, 10, 15, 4, 26, 68, 92, 23, 60, 83, 59, 82, 38,
                       45, 74, 48, 20, 3, 2, 13, 34, 57, 87, 100, 32, 30, 56, 70, 75, 44, 12,
                       93, 86, 33, 5, 53, 49, 71, 42, 6, 41, 16, 80, 78, 39, 24, 76, 8, 69, 47,
                       55, 19, 25, 37, 99, 67, 27, 51, 9, 14, 91, 54, 79, 81, 73, 66, 22, 85, 11,
                       1, 29, 96, 35, 98, 50, 18, 95, 90, 52, 40, 17, 72, 63, 21, 28, 64)


getTaxaColumnNum <- function(myT)
{
	colNames <- names(myT)

	for( i in 1:length(colNames))
	{
		if( grepl( "Imagination", colNames[i]))
			return (i + 1);
	}

	return (-1);
}

for ( t in taxa )
{
	inFileName <- paste( "pivoted_", t, "asColumnsLogNormalPlusMetadata.txt" , sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",2),"numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
	myT <- myT[ myT$read== "r1" & ! is.na(myT$AGE) , ]
        myT <- myT[as.numeric(sub("r1_", "", myT$id)) %in% selectedSamplesPC,]


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

	pdf( paste(t , "highabundance_samp_cohort_correlationPlots.pdf",sep=""))

	taxCol <- getTaxaColumnNum(myT)

        myT$cohortNum <- as.numeric(myT$Study.ID <= "HC47")

      	inFileName <- paste(t, "_SparseThreeCol.txt", sep = "")
	myCountT <- read.table(inFileName, header=FALSE, sep = "\t")
	myColClasses <- c(rep("character", 2), "numeric")
	myCountT <- read.table(inFileName, header = FALSE, sep = "\t", colClasses = myColClasses)
        colnames(myCountT) <- c("id", "taxa", "count")
        myCountT$id <- unlist(lapply(strsplit(myCountT$id, split="\\."), "[[", 1))
        myCountT$read <- unlist(lapply(strsplit(myCountT$id, split="_"), "[[", 1))
	myCountT <- myCountT[ myCountT$read == "r1", ]

        highAbundance<-vector()
        for(i in myCountT$id){
            highAbundance[i] <- sum(myCountT$count[which(myCountT$id == i, arr.ind=TRUE)])
        }
        highAbundance<-as.data.frame(highAbundance)
        highAbundance$id <- rownames(highAbundance)
        merged_df <- merge(myCountT, highAbundance)
        merged_df$percent <- merged_df$count / merged_df$highAbundance
        keptatnPercent <- unique(merged_df$taxa[which(merged_df$percent > nPercent)])

        #Reduce the number of columns by one to account for the cohort column
	for( i in c(3, taxCol : (ncol(myT) - 1) ))
	{
		if( sum( myT[,i] >0 ) > nrow(myT) /4 && names(myT)[i] %in% keptatnPercent)
                    {

			 for ( j in 5:(taxCol-1))
			 {
			 	namesA[index] <- names(myT)[i]
			 	namesB[index] <- names(myT)[j]

			 	rVals[index] <- cor( myT[,i], myT[,j], method="spearman")
			 	aLm <- lm(myT[,i] ~ myT[,j])
			 	pValueLinear[index] <- anova(aLm)$"Pr(>F)"[1]
			 	kendallP[index] <- Kendall(myT[,i], myT[,j])$sl[1]
                                cohortNum <- myT$cohortNum

                                bug <- myT[,i]
                                testvar <- myT[,j]

                                myFrame <- data.frame(bug, testvar, cohortNum)

                                fullModel <- gls(bug ~ testvar, method="REML", correlation=corCompSymm(form=~1|factor(cohortNum)), data = myFrame )
                                reducedModel <- gls( bug ~ testvar, method="REML", data = myFrame)
                                fullModelLME <- lme( bug ~ testvar, method="REML", random = ~1|factor(cohortNum), data = myFrame)

                                pValueTestvar[index] <- anova(fullModelLME)$"p-value"[2]
                                pValueCohort[index] <- anova(fullModelLME, reducedModel)$"p-value"[2]
                                iccCohort[index] <- coef(fullModel$modelStruct[1]$corStruct, unconstrained=FALSE)[[1]]

			 	myText <- paste( namesA[index] ,namesB[index] ,"\n", "p=", format(pValueLinear[index] , digits=3),
                                                "r=", format( rVals[index], digits=3),
                                                "kendall p=" ,
                                                format( kendallP[index], digits=3),
                                                "\n p_Testvar=", format(pValueTestvar[index], digits=3),
                                                "p_cohort=", format(pValueCohort[index], digits=3),
                                                "icc Cohort=", format(iccCohort[index], digits=3)
                                                )
			 		plot(myT[,j],myT[,i] , main=myText, ylab =namesA[index], xlab=namesB[index]  )
			 		abline(aLm)

			 	index <- index + 1
			 }
		}
	}

dev.off()

dFrame <- data.frame( namesA, namesB, kendallP, pValueLinear, rVals, pValueTestvar, pValueCohort, iccCohort)

dFrame <- dFrame [order(dFrame$kendallP),]
dFrame$adjKendall <-  p.adjust( dFrame$kendallP, method = "BH" )
dFrame$adjPLinear <-  p.adjust( dFrame$pValueLinear, method = "BH" )
dFrame$adjTestvar <- p.adjust( dFrame$pValueTestvar, method = "BH" )
dFrame$adjCohort <- p.adjust(dFrame$pValueCohort, method = "BH" )

write.table( file= paste("highabundance_samp_cohort_pValuesTaxaVsMetadata_", t, "_", nPercent, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}

