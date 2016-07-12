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
	names <- vector()
	pValuesPatientID <- vector()
	pValuesCalorimetry <- vector()
        pValuesTime <- vector()
        pValuesBug <- vector()
        pValuesIntake <- vector()
        pValuesInteraction <- vector()
	meanBug <- vector()

	index <- 1

        pdf( paste(taxa, "_BothTime_Relative_Function.pdf", sep=""))

	for( i in 3:(numCols - 7))
		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			#bug <- log10( myT[,i] + 0.00001)
			#meanBug[index] <- mean(bug)
                        bug <- myT[,i]
                        time <- factor(myT$Time)
			patientID <- myT$Sample
			## calorimetry<- myT$Relative.Energy.Content
                        calorimetry <- myT$cal.g
                        intake <- myT$Energy.Intake

			myFrame <- data.frame(bug, time, patientID, calorimetry, intake)

			fullModel <- lm( calorimetry~  bug*intake + time)
                        			fullModel <- gls( calorimetry~  bug*intake + intake,
				 method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
				data = myFrame )

			reducedModel <- gls( calorimetry~  bug*intake + time, method="REML",	data = myFrame )

			fullModelLME <- lme( calorimetry~  bug*intake + time, method="REML", random = ~1|factor(patientID), data = myFrame)

			pValuesBug[index] <- anova(fullModelLME)$"p-value"[2]
			pValuesTime[index] <- anova(fullModelLME)$"p-value"[4]
                        pValuesIntake[index] <- anova(fullModelLME)$"p-value"[3]
                        pValuesInteraction[index] <- anova(fullModelLME)$"p-value"[5]
			pValuesPatientID[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

			## pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], " pTime=", format(pValuesTime[index], digits=3), "\n",
                            " pValuesBug= ", format(pValuesBug[index],digits=3), "\n",
                            " pValuesIntake= ", format(pValuesIntake[index], digits=3), "\n",
                            " pValuesInteraction= ", format(pValuesInteraction[index], digits=3),
									" pValuesPatientID= " , format(	pValuesPatientID[index], digits=3), "\n",
										" icc= " , format( intraclassCoefficient, digits=3 ), sep="")

			plot( calorimetry ~ bug + time*intake, ylab = names[index],
					main = graphMain )

			index=index+1

		}
	dFrame <- data.frame( names, pValuesBug, pValuesTime, pValuesIntake, pValuesInteraction, pValuesPatientID, intraclassCoefficient)
	dFrame <- dFrame [order(dFrame$pValuesIntake),]
        dFrame$adjustedpValuesBug <- p.adjust( dFrame$pValuesBug, method = "BH" )
        dFrame$adjustedpValuesTime <- p.adjust( dFrame$pValuesTime, method = "BH" )
        dFrame$adjustedpValuesIntake <- p.adjust( dFrame$pValuesIntake, method = "BH" )
        dFrame$adjustedpValuesInteraction <- p.adjust( dFrame$pValuesInteraction, method = "BH")
        dFrame$adjustedpValuesPatientID <- p.adjust( dFrame$pValuesPatientID, method = "BH" )

	write.table(dFrame, file=paste("pValuesFor_", taxa, "_CalorimetryVS_BugIntake_Time_Patient_RELATIVE_function.txt",sep=""), sep="\t",row.names=FALSE)
		dev.off()
}
