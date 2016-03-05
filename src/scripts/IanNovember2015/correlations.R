rm(list=ls())

# setwd("C:\\Caroll_Nov_2015\\spreadsheets")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Anorexia")

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("vegan")

taxa <- c("phylum","class","order","family","genus")

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

	pdf( paste(t , "NEWcorrelationPlots.pdf",sep=""))

	taxCol <- getTaxaColumnNum(myT)

        myT$cohortNum <- as.numeric(myT$Study.ID <= "HC50")

        #Reduce the number of columns by one to account for the cohort column
	for( i in c(3, taxCol : (ncol(myT) - 1) ))
	{
		if( sum( myT[,i] >0 ) > nrow(myT) /4 )
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

write.table( file= paste( "NEWpValuesTaxaVsMetadata_", t, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}

