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
                                        #t<-taxaLevels[[1]]
        pdf( paste(t, "_DayPatient_NoLow_plots.pdf", sep = ""))
      	inFileName <- paste(t, "LogNormalwithMetadataDailyR2_Edit.txt", sep="")
	myT <-read.table(inFileName, header=TRUE, sep="\t")
	numCols <- ncol(myT)
        # I don't think that will necessarily be right.
	myColClasses <- c(rep("character", 2), "numeric", "character", rep("numeric", numCols-4))
	myT <-read.table(inFileName, header=TRUE, sep="\t", colClasses=myColClasses)
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
                        EnergyIntake<-myT$Energy.Intake..kcal.day.
                        taxaType <- as.numeric(myT[,i])

                        DayPatient<-lm(taxaType ~  Day*patient, x = TRUE)

                        DayPatientpVal[index] <- list(summary(DayPatient)$coefficients[,4][-1])

                        #Compute Shannon diversity and Shannon richness via vegan

                        #myShannon <- diversity(myT[,i])
                        #myRichness <-
                        # Second order in time
                        #time2ndOrderPseud <- lm( taxaType ~ poly( Day, 2), x=TRUE)
                        #This is probably just developed for the simple regular intervals of time points.

                        # Autocorrelation Code
                        #M1 <- gls( taxaType ~ EnergyIntake , data= myFrame )
                        #M2 <- gls( taxaType ~ EnergyIntake + Day , data= myFrame )
                        #M3 <- gls( taxaType ~ EnergyIntake + Day , data= myFrame , correlation = corCompSymm(form=~Day ))
                        #M4 <- gls( taxaType ~ EnergyIntake + Day , data= myFrame , correlation = corAR1(form=~Day ))
                        #M5 <- gls( taxaType ~ EnergyIntake * Day , data= myFrame  )

                        #M6 <- gls( taxaType ~ EnergyIntake + poly(Day,2) , data= myFrame )
                        #M7 <- gls( taxaType ~ EnergyIntake + poly(Day,2) , data= myFrame ,correlation = corCompSymm(form=~Day ))
                        #M8 <- gls( taxaType ~ EnergyIntake + poly(Day,2) , data= myFrame ,correlation = corAR1(form=~Day ))

                        #E <- residuals( M3, type="normalized")
                        #acf(E)
                        #anova(M2,M3,M4, M5,M6)
                        #anova(M6,M7,M8)
                        #summary(M2)
                        #summary(M4)
                        #plot(residuals(M2), myFrame$Day)
                        #boxplot( residuals(M2) ~ myFrame$EnergyIntake)

                        # Compiling the p-values for eventual print out.

                        names[index] = names(myT)[i]

                                        #Graphs for each of the models here...
                        graphMain = paste(names(myT)[i], "\n",
                            "pDay=", format(DayPatientpVal[[index]][1], digits=3), "\n",
                        "pPatientB=", format(DayPatientpVal[[index]][2], digits=3),
                            "pPatientC=", format(DayPatientpVal[[index]][3], digits=3), "\n",
                            "pEnergyIntake=", format(DayPatientpVal[[index]][4], digits=3), "\n",
        "pDay:PatientB=", format(DayPatientpVal[[index]][5], digits=3),
"pDay:PatientC=", format(DayPatientpVal[[index]][6], digits=3))
par(mar = c(5, 4, 6, 2))
                        plot(EnergyIntake, taxaType, col=colors, main=graphMain)
                        abline(a = DayPatient$coef[1], b = DayPatient$coef[2])
                        abline(a = DayPatient$coef[1] + DayPatient$coef[3], b = DayPatient$coef[5], col="BLUE")
                        abline(a = DayPatient$coef[1] + DayPatient$coef[4], b = DayPatient$coef[6], col="RED")
                       index = index + 1
		}
            }


        # Building the data.frames to eventually print out the p-values
        DayPatientPV.df <- data.frame(DayPatientpVal)
        DayPatientPV.df <- t(DayPatientPV.df)
        dFrameDayPatient <- data.frame(names, DayPatientPV.df)

        for (m in 2:dim(dFrameDayPatient)[2])
        {
           dFrameDayPatient[,dim(dFrameDayPatient)[2] + 1] <- p.adjust(dFrameDayPatient[,m], method = "BH")
           colnames(dFrameDayPatient)[ncol(dFrameDayPatient)]<-paste0("adj",colnames(dFrameDayPatient)[m])
        }

                                        #Finally, writing out the p-values and BH adjusted p-values
        write.table(dFrameDayPatient, file = paste("pValuesLongPatient_Day_NoLow_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

    dev.off()
}
