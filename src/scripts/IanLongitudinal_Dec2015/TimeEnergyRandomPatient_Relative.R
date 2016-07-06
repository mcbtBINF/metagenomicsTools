rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")
library("Kendall")
library("vegan")
library("lmtest")
library("pscl")
library("nlme")
library("gtools")

taxaLevels <- c("phylum","class","order","family","genus")

for(t in taxaLevels )
{
      	inFileName <- paste(t,"_LogNormalwithMetadataWeekly_NearestSamplewithDay_withCals.txt", sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",2),"numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
	myT <- myT[ !is.na(myT[2]), ]

	patientPValues <- vector()
	timePValues <- vector()
	interactionPValues <- vector()
	names <-vector()

        BMIPValues <- vector()
        REEPValues <- vector()
        DITPValues <- vector()
        PAEPValues <- vector()
        PAENormPValues <- vector()
        AEEPValues <- vector()
        energyPValues <- vector()
        patient <- vector()
        EnergyContentPValues <- vector()

        myLmpVal <- list()

	colors <- vector()
	cIndex <-1
	for ( j in 1: nrow(myT))
	{
		if( substr(myT[j,]$Sample.ID,1,1) == "B")
		{
                    colors[cIndex] <- "Blue"
                    patient[cIndex] <- "B"
		} else if (substr(myT[j,]$Sample.ID,1,1) == "C" )
		{
                    colors[cIndex] <- "Red"
                    patient[cIndex] <- "C"
		} else
		{
                    colors[cIndex] <- "Black"
                    patient[cIndex] <- "A"
		}
		cIndex = cIndex + 1

	}

	index <-1
        selectRows <- which(is.na(myT$Day) == FALSE)
        myT <- myT[selectRows,]
        colors <- colors[selectRows]
        patient <- patient[selectRows]
        myT <- myT[mixedorder(myT[,1]),]

        ## myT$Shannon <- apply(myT[,2:(ncol(myT)-15)],1,diversity)
        ## Shannon <- myT$Shannon
        ## Day <- myT$Day
        ## EnergyContent <- myT$Energy.Content..cal.g/myT$Energy.Intake..kcal.day.

##        ImputedBMI<-myT$Imputed.BMI
##        EnergyIntake<-myT$Energy.Intake..kcal.day.


        ## Shannonlm <- lm(Shannon ~ patient*EnergyContent, x = TRUE)
        ## #ShannonlmBMI <- lm(Shannon ~ ImputedBMI*patient, x = TRUE)
        ## #ShannonlmEnergy <- lm(Shannon ~ EnergyIntake*patient, x = TRUE)
        ## ShannonP[indexS] <- list(summary(Shannonlm)$coefficients[,4][-1])

        ## indexS <- indexS + 1

        ## # Building the data.frames to eventually print out the p-values
        ## DayPatientPV.df <- data.frame(ShannonP)
        ## DayPatientPV.df <- t(DayPatientPV.df)
        ## dFrameDayPatient <- DayPatientPV.df

        ## plot(Day, Shannon, col=colors)

        ## write.table(dFrameDayPatient, file = paste("pValuesShannon_DayPatient_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

        ## dev.off()

	pdf(paste(t, "_TimeEnergyRandomPatient_RELATIVE_plots.pdf", sep="") )

        for( i in 2:(ncol(myT) - 15))
            {
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                        bug <- myT[,i]
                        calorimetry <- myT$Energy.Content..cal.g/myT$Energy.Intake..kcal.day.
                        time <- myT$Day
                        ## time <- myT$Day %/% 7
                        patientID <- as.numeric(factor(patient))

			myFrame <- data.frame(bug, time, patientID, calorimetry)
                        names[index] <- names(myT)[i]

       			fullModel <- gls( bug~  calorimetry + time,
		 method="REML", correlation=corCompSymm(form=~1|factor(patientID)),
				data = myFrame )

			reducedModel <- gls( bug~  calorimetry + time, method="REML",	data = myFrame )

			fullModelLME <- lme(bug~  calorimetry + time, method="REML", random = ~1|factor(patientID), data = myFrame)

			pValuesCalorimetry[index] <- anova(fullModelLME)$"p-value"[2]
			pValuesTime[index] <- anova(fullModelLME)$"p-value"[3]
			pValuesPatientID[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

                        meanA <- mean(as.numeric(myT[colors == "Black", i]))
                        meanB <- mean(as.numeric(myT[colors == "Blue", i]))
                        meanC <- mean(as.numeric(myT[colors == "Red", i]))

			graphMain =  paste( names(myT)[i], " pTime=", format(pValuesTime[index], digits=3), "\n",
								" pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
									" pValuesPatientID= " , format(	pValuesPatientID[index], digits=3), "\n",
										" icc= " , format( intraclassCoefficient, digits=3 ), sep="")


                        plot(time, bug, main = graphMain, col=colors, xlab="Time (Days)", ylab = "Log-Abundance", pch=16, ylim=c(0,5))

                        abline(a = meanA
                              #     myLm$coef[1]
                             , b = myLm$coef[4])
                        abline(a = meanB
                                        #myLm$coef[1] + myLm$coef[2] #+ myLm$coef[6]
                             , b = myLm$coef[4] + myLm$coef[5] #+ myLm$coef[2]
                             , col="BLUE")
                        abline(a = meanC
                                   #myLm$coef[1] + myLm$coef[3] #+ myLm$coef[4]
                             , b = myLm$coef[4] + myLm$coef[6] #+ myLm$coef[2]
                             , col="RED")

                     index = index + 1
		}
	}
	dev.off()

	dFrame <- data.frame( names, pValuesCalorimetry, pValuesTime, pValuesPatientID, intraclassCoefficient)
	dFrame <- dFrame [order(dFrame$pValuesCalorimetry),]
	dFrame$adjustedpValuesCalorimetry <- p.adjust( dFrame$pValuesCalorimetry, method = "BH" )
        dFrame$adjustedpValuesTime <- p.adjust( dFrame$pValuesTime, method = "BH" )
        dFrame$adjustedpValuesPatientID <- p.adjust( dFrame$pValuesPatientID, method = "BH" )

	write.table(dFrame, file= paste( "TimeEnergyRandomPatient_RELATIVE_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

#        dev.off()
}
