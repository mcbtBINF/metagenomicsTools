rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

## setwd("C:\\Susan_Oct2015")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/BombCalorimetry")
taxaLevels <- c("phylum","class","order","family","genus")

for(taxa in taxaLevels )
{
	inFileName <- paste(taxa, "_paired_metadata.txt", sep ="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	## myColClasses <- c(rep("character",1), rep("numeric", numCols-1))
        myColClasses <- c(rep("numeric", numCols - 6), "character", rep("numeric", 5))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

	myT <- myT[ ! is.na(myT$cal.g), ]
        ## This drops sample 712 from consideration.

        pairedList <- c(701, 702, 703, 705, 706, 707, 708, 709, 710)
        ## to paired samples only
        myT <- myT[myT$Sample %in% pairedList,]

	names <- vector()
	pValuesTime <- vector()
	pValuesSubject <- vector()
	pValuesCalorimetry <- vector()
        pValuesBMI <- vector()
        pValuesLOS <- vector()
        pValuesAge <- vector()
        pValuesEnergy.Intake <- vector()
        pValuesRelative.Energy.Content <- vector()

	meanBug <- vector()
	index <- 1
	pdf( paste(taxa, "_new_plots.pdf", sep=""))

	for( i in 3:(numCols - 7))
		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
                    {
                        ## Why is this log10'd?
			## bug <- log10( myT[,i] + 0.00001)
                        bug <- myT[,i]
                        meanBug[index] <- mean(bug)
			time <- factor(myT$Time)
			patientID <- myT$Sample
			calorimetry<- myT$cal.g
                        age <- myT$Age
                        BMI <- myT$BMI
                        race <- myT$Race
                        LOS <- myT$LOS
                        EnergyIntake <- myT$Energy.Intake
                        RelativeEnergyContent <- myT$Relative.Energy.Content

			myFrame <- data.frame(bug, time, patientID, calorimetry, age, BMI, race, LOS, EnergyIntake, RelativeEnergyContent)

			fullModel <- gls( bug~  time + calorimetry,
				 method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
				data = myFrame )

			reducedModel <- gls( bug~  time + calorimetry, method="REML",	data = myFrame )

			fullModelLME <- lme(bug~  time + calorimetry, method="REML", random = ~1|factor(patientID), data = myFrame)

			pValuesTime[index] <- anova(fullModelLME)$"p-value"[2]
			pValuesCalorimetry[index] <- anova(fullModelLME)$"p-value"[3]
			pValuesSubject[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], " pTime=", format(pValuesTime[index], digits=3), "\n",
								" pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
									" pSubject= " , format(	pValuesSubject[index], digits=3), "\n",
										" icc= " , format( intraclassCoefficient, digits=3 ), sep="")

			plot( bug ~ calorimetry, ylab = names[index],
					main = graphMain )
			index=index+1

		}

	dFrame <- data.frame( names, pValuesTime ,pValuesSubject,pValuesCalorimetry ,meanBug)
	dFrame <- dFrame [order(dFrame$pValuesCalorimetry),]
        dFrame$adjTime <- p.adjust(dFrame$pValuesTime, method="BH")
        dFrame$adjSubject <- p.adjust(dFrame$pValuesSubject, method="BH")
	dFrame$adjustedpValuesCalorimetry <- p.adjust( dFrame$pValuesCalorimetry, method = "BH" )
	write.table(dFrame, file=paste("NEW_pValuesFor", taxa, ".txt",sep=""), sep="\t",row.names=FALSE)
		dev.off()
}
