#Using BMI rather than Imputed.BMI

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

for(t in taxaLevels )
{
        pdf( paste(t, "_BMIPatient_ANOVA_NoLow_plots.pdf", sep = ""))
      	inFileName <- paste(t, "LogNormalwithMetadataDailyR2_Edit.txt", sep="")
	myT <-read.table(inFileName, header=TRUE, sep="\t")
	numCols <- ncol(myT)
        # I don't think that will necessarily be right.
	myColClasses <- c(rep("character", 2), "numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName, header=TRUE, sep="\t", colClasses=myColClasses)
	myT <- myT[ !is.na(myT[2]), ]

	names <-vector()

        BMIPatientpVal <- list()
        OLDBMIPatientpVal <- list()

        #Patients to colors
	colors <- vector()
        patient <- vector()
	cIndex <- 1

        myT <- myT[-which(myT$Sample.ID %in% list(37, 45, 52, 58, 7, 9),arr.ind=TRUE),]

        # Dropping so that we only consider the cases with BMI present (conservative)
        myT <- myT[which(is.na(myT$BMI) == FALSE, arr.ind=TRUE),]

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

        # Should reliably do this even for the mixed case
        myT <- myT[mixedorder(myT[,1]),]
        #Remove the low depth samples.
        for( i in 2:(ncol(myT) - 14) )
            {
                #Remove from consideration rare organisms.
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                        #Makes these components easier to work with and compute derivative values
                        Day<-myT$Day
                        ImputedBMI<-myT$Imputed.BMI
                        BMI <- myT$BMI
                        EnergyIntake<-myT$Energy.Intake..kcal.day.
                        taxaType <- as.numeric(myT[,i])

                        BMIPatient<-lm(taxaType ~  BMI*patient, x = TRUE)

                        OLDBMIPatientpVal[index] <- list(summary(BMIPatient)$coefficients[,4][-1])
                        BMIPatientpVal[index] <- list(anova(BMIPatient)$"Pr(>F)"[1:3])

                        # Compiling the p-values for eventual print out.

                        names[index] = names(myT)[i]

                        #Graphs for each of the models here...
                        # These are uncorrected p-values
                        graphMain = paste(names(myT)[i]), "\n",
#                            "pBMI=", format(OLDBMIPatientpVal[[index]][1], digits=3), "\n",
#                        "pPatientB=", format(OLDBMIPatientpVal[[index]][2], digits=3),
#                            "pPatientC=", format(OLDBMIPatientpVal[[index]][3], digits=3), "\n",
#                           # "pEnergyIntake=", format(OLDBMIPatientpVal[[index]][4], digits=3), "\n",
#        "pBMI:PatientB=", format(OLDBMIPatientpVal[[index]][4], digits=3),
#"pBMI:PatientC=", format(OLDBMIPatientpVal[[index]][5], digits=3))
par(mar = c(5, 4, 6, 2))
                        plot(BMI, taxaType, col=colors, main=graphMain)
                        abline(a = BMIPatient$coef[1], b = BMIPatient$coef[2])
                        abline(a = BMIPatient$coef[1] + BMIPatient$coef[3], b = BMIPatient$coef[5] + BMIPatient$coef[2], col="BLUE")
                        abline(a = BMIPatient$coef[1] + BMIPatient$coef[4], b = BMIPatient$coef[6] + BMIPatient$coef[2], col="RED")
                       index = index + 1
		}
            }

        # Building the data.frames to eventually print out the p-values
        BMIPatientPV.df <- data.frame(BMIPatientpVal)
        BMIPatientPV.df <- t(BMIPatientPV.df)
        modeldf <- as.data.frame(matrix(unlist(OLDBMIPatientpVal), nrow=length(OLDBMIPatientpVal), byrow = TRUE))
        dFrameBMIPatient <- data.frame(names, BMIPatientPV.df, modeldf)
        colnames(dFrameBMIPatient) <- c("names", "ANOVA->BMI", "ANOVA->patient", "ANOVA->BMI:patient", "BMI", "patientB", "patientC", "BMI:patientB", "BMI:patientC")
        for (m in 2:dim(dFrameBMIPatient)[2])
        {
           dFrameBMIPatient[,dim(dFrameBMIPatient)[2] + 1] <- p.adjust(dFrameBMIPatient[,m], method = "BH")
           colnames(dFrameBMIPatient)[ncol(dFrameBMIPatient)]<-paste0("adj",colnames(dFrameBMIPatient)[m])
        }

        #Finally, writing out the p-values and BH adjusted p-values
        write.table(dFrameBMIPatient, file = paste("pValuesLongPatient_BMI_ANOVA_NoLow_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

    dev.off()
}
