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
        pdf( paste(t, "_FinalR1_DayPatient_ANOVA_NoLow_plots.pdf", sep = ""))
      	inFileName <- paste(t, "LogNormalwithMetadata_Edit.txt", sep="")
        ## It was "LogNormalwithMetadataDailyR2_Edit.txt"
	myT <-read.table(inFileName, header=TRUE, sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character", 2), "numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName, header=TRUE, sep="\t", colClasses=myColClasses)
	myT <- myT[ !is.na(myT[2]), ]

	names <-vector()
        ## For swapping and use in the plotting of corrected p-values.
        DayPatientpVal <- list()
        OLDDayPatientpVal <- list()

        ## Patients to colors
	colors <- vector()
        patient <- vector()
	cIndex <- 1
        ## Removes samples of low sequencing depth.
        ## Will need to check on this again.
        myT <- myT[-which(myT$Sample.ID %in% list(37, 45, 52, 58, 7, 9),arr.ind=TRUE),]

        ## Create a color vector for graphical output according to patient letter.
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

        # Reorder data
        myT <- myT[mixedorder(myT[,1]),]
        ## 14 is the offset due to the metadata columns
        for( i in 2:(ncol(myT) - 14) )
            {
        ## Remove from consideration rare organisms.
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                        ## Easy to refer to variable names
                        Day<-myT$Day
                        ImputedBMI<-myT$Imputed.BMI
                        BMI <- myT$BMI
                        EnergyIntake<-myT$Energy.Intake..kcal.day.
                        taxaType <- as.numeric(myT[,i])

                        ## Modeling

                        DayPatient<-lm(taxaType ~  Day*patient, x = TRUE)

                        ##Store coefficients for plotting lines

                        OLDDayPatientpVal[index] <- list(summary(DayPatient)$coefficients[,4][-1])
                        ## Store p-values for subsequent correction
                        DayPatientpVal[index] <- list(anova(DayPatient)$"Pr(>F)"[1:3])

                        names[index] = names(myT)[i]
                        index = index + 1
                    }
            }

        # Building the data.frames to eventually print out the p-values
        DayPatientPV.df <- data.frame(DayPatientpVal)
        DayPatientPV.df <- t(DayPatientPV.df)
        modeldf <- as.data.frame(matrix(unlist(OLDDayPatientpVal), nrow=length(OLDDayPatientpVal), byrow = TRUE))
        dFrameDayPatient <- data.frame(names, DayPatientPV.df, modeldf)
        colnames(dFrameDayPatient) <- c("names", "ANOVA->Day", "ANOVA->patient", "ANOVA->Day:patient", "Day", "patientB", "patientC", "Day:patientB", "Day:patientC")

        for (m in 2:dim(dFrameDayPatient)[2])
        {
           dFrameDayPatient[,dim(dFrameDayPatient)[2] + 1] <- p.adjust(dFrameDayPatient[,m], method = "BH")
           colnames(dFrameDayPatient)[ncol(dFrameDayPatient)]<-paste0("adj",colnames(dFrameDayPatient)[m])
        }

        # Repeat modeling so as to use corrected p-values for the graph display
	index <-1

        myT <- myT[mixedorder(myT[,1]),]
        for( i in 2:(ncol(myT) - 14) )
            {
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                        Day<-myT$Day
                        ImputedBMI<-myT$Imputed.BMI
                        BMI <- myT$BMI
                        EnergyIntake<-myT$Energy.Intake..kcal.day.
                        taxaType <- as.numeric(myT[,i])

                        DayPatient<-lm(taxaType ~ Day*patient, x = TRUE)

                        OLDDayPatientpVal[index] <- list(summary(DayPatient)$coefficients[,4][-1])
                        DayPatientpVal[index] <- list(anova(DayPatient)$"Pr(>F)"[1:3])

                        names[index] = names(myT)[i]

                        ## Graphs for each of the models here...
                        ## These are corrected p-values now
                        graphMain = paste(names(myT)[i]) ## , "\n",
                        ## "pDay=", format(dFrameDayPatient[index, 13], digits=3), "\n",
                        ## "pPatientB=", format(dFrameDayPatient[index, 14], digits=3),
                        ## "pPatientC=", format(dFrameDayPatient[index, 15], digits=3), "\n",
                        ## "pEnergyIntake=", format(dFrameDayPatient(index,[4], digits=3), "\n",
                        ## "pDay:PatientB=", format(dFrameDayPatient[index,16], digits=3),
                        ## "pDay:PatientC=", format(dFrameDayPatient[index, 17], digits=3))
                        par(mar = c(5, 4, 6, 2))

                        plot(Day, taxaType, col=colors, main=graphMain, pch=16, xlab="Days Under Treatment", ylab = "Log Normalized Abundance")
                        legend("bottomright",
                               c("Patient A", "Patient B", "Patient C"),
                               pch = c(16, 16, 16),
                               col=c("BLACK", "BLUE", "RED"))
                        ## abline(a = DayPatient$coef[1], b = DayPatient$coef[2])
                        ## abline(a = DayPatient$coef[1] + DayPatient$coef[3], b = DayPatient$coef[5] + DayPatient$coef[2], col="BLUE")
                        ## abline(a = DayPatient$coef[1] + DayPatient$coef[4], b = DayPatient$coef[6] + DayPatient$coef[2], col="RED")

                        index = index + 1
		}
            }

        #Finally, writing out the p-values and BH adjusted p-values
        write.table(dFrameDayPatient, file = paste("FinalR1_pValuesLongPatient_Day_ANOVA_NoLow_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

    dev.off()
}
