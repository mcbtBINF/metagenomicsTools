rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")

library("Kendall")
library("vegan")
library("lmtest")
library("pscl")
library("nlme")
library("gtools")

taxaLevels <- c("phylum","class","order","family","genus")
#Right Now I'm being lazy with the changing of the variable names
indexS <- 1
ShannonP <- list()
ShannonSummary <- list()
for(t in taxaLevels )
{
                                        #t<-taxaLevels[[1]]
        pdf( paste("Shannon_","ImputedBMIPatient_",t, "_plots.pdf", sep = ""))
      	inFileName <- paste(t, "LogNormalwithMetadataDailyR2_Edit.txt", sep="")
	myT <-read.table(inFileName, header=TRUE, sep="\t")
	numCols <- ncol(myT)
        # I don't think that will necessarily be right.
	#myColClasses <- c(rep("character", 2), "numeric", "character", rep("numeric", numCols-4))
	#myT <-read.table(inFileName, header=TRUE, sep="\t", colClasses=myColClasses)
	myT <- myT[ !is.na(myT[2]), ]

	names <-vector()

        DayPatientpVal <- list()

        #Patients to colors
	colors <- vector()
        patient <- vector()
	cIndex <- 1

         myT <- myT[-which(myT$Sample.ID %in% list(37, 45, 52, 58, 7, 9),arr.ind=TRUE),]

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

#	index <-1

        # Should reliably do this even for the mixed case
        myT <- myT[mixedorder(myT[,1]),]
        myT$Shannon <- apply(myT[,2:(ncol(myT)-14)],1,diversity)
        Shannon <- myT$Shannon
        Day <- myT$Day
        ImputedBMI<-myT$Imputed.BMI
        EnergyIntake<-myT$Energy.Intake..kcal.day.

        Shannonlm <- lm(Shannon ~ ImputedBMI*patient, x = TRUE)
        #ShannonlmBMI <- lm(Shannon ~ ImputedBMI*patient, x = TRUE)
        #ShannonlmEnergy <- lm(Shannon ~ EnergyIntake*patient, x = TRUE)
        ShannonP[indexS] <- list(summary(Shannonlm)$coefficients[,4][-1])
#        ShannonSummary[indexS]<-summary(Shannonlm)

#        names[indexS] = names(myT)[i]

        indexS <- indexS + 1

        # Building the data.frames to eventually print out the p-values
        DayPatientPV.df <- data.frame(ShannonP)
        DayPatientPV.df <- t(DayPatientPV.df)
        dFrameDayPatient <- DayPatientPV.df

        plot(Day, Shannon, col=colors)

        write.table(dFrameDayPatient, file = paste("pValuesShannon_ImputedBMIPatient_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

        dev.off()

#        dFrameDayPatient <- data.frame(names, DayPatientPV.df)

#        for (m in 2:dim(dFrameDayPatient)[2])
#        {
#           dFrameDayPatient[,dim(dFrameDayPatient)[2] + 1] <- p.adjust(dFrameDayPatient[,m], method = "BH")
           #colnames(dFrameDayPatient)[ncol(dFrameDayPatient)]<-paste0("adj",colnames(dFrameDayPatient)[m])
#        }

                                        #Finally, writing out the p-values and BH adjusted p-values
#        write.table(dFrameDayPatient, file = paste("pValuesLongDayPatientShannon_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

#    dev.off()
}
