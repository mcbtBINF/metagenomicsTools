rm(list=ls())
library("lmtest")
library("nlme")
library("pscl")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/BombCalorimetry/")

taxaLevels <- c("phylum","class","order","family","genus")

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa,  "_paired_metadata.txt", sep ="")
	myT <-read.table(inFileName, header=TRUE,sep="\t")
	numCols <- ncol(myT)
##	myColClasses <- c(rep("character",1), rep("numeric", numCols-1))
##	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

##	myT <- myT[ myT$timepoint == "2" &  ! is.na(myT$calorimetryData), ]
        myT <- myT[myT$Sample != 712,]
        myT <- myT[myT$Time == 2,]
	names <- vector()
	pValuesPatientID <- vector()
	pValuesCalorimetry <- vector()
        pValuesTime <- vector()
	meanBug <- vector()

	index <- 1

        pdf( paste(taxa, "_Time2_Relative.pdf", sep=""))

	for( i in 3:(numCols - 7))
		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			#bug <- log10( myT[,i] + 0.00001)
			#meanBug[index] <- mean(bug)
                        bug <- myT[,i]
                        time <- factor(myT$Time)
			patientID <- myT$Sample
			calorimetry<- myT$Relative.Energy.Content

			myFrame <- data.frame(bug, time, patientID, calorimetry)

			fullModel <- lm( bug~  calorimetry + patientID)
#      			fullModel <- gls( bug~  calorimetry,				 method="REML",correlation=corCompSymm(form=~1|factor(patientID)), 				data = myFrame )

#			reducedModel <- gls( bug~  calorimetry, method="REML",	data = myFrame )

#			fullModelLME <- lme(bug~  calorimetry, method="REML", random = ~1|factor(patientID), data = myFrame)

			pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
                        pValuesPatientID[index] <- anova(fullModel)$"Pr(>F)"[2]
                                        #			pValuesTime[index] <- anova(fullModelLME)$"p-value"[3]

                                        #			pValuesPatientID[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
##			intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


			## pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], "\n",
								" pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
									" pValuesPatientID= " , format(	pValuesPatientID[index], digits=3))#, ##"\n",
##										" icc= " , format( intraclassCoefficient, digits=3 ), sep="")
#x)
			plot( bug ~ calorimetry, ylab = names[index],
					main = graphMain )

			index=index+1

		}
	dFrame <- data.frame( names, pValuesCalorimetry, pValuesPatientID) ## ,intraclassCoefficient)
	dFrame <- dFrame [order(dFrame$pValuesCalorimetry),]
	dFrame$adjustedpValuesCalorimetry <- p.adjust( dFrame$pValuesCalorimetry, method = "BH" )
#        dFrame$adjustedpValuesTime <- p.adjust( dFrame$pValuesTime, method = "BH" )
        dFrame$adjustedpValuesPatientID <- p.adjust( dFrame$pValuesPatientID, method = "BH" )

	write.table(dFrame, file=paste("pValuesFor_", taxa, "_Calorimetry_Patient_Time2_RELATIVE.txt",sep=""), sep="\t",row.names=FALSE)
		dev.off()
}
