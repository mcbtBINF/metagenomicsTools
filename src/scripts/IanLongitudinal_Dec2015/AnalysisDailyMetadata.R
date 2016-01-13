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
        pdf( paste(t, "_plots.pdf", sep = ""))
      	inFileName <- paste(t,"LogNormalwithMetadataDailyR2_Edit.txt", sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
        # I don't think that will necessarily be right.
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

        fullPValues <- list()
        noIntPValues <- list()
        BMIPValues <- list()
        EnergyPValues <- list()
        DayBMIPV <- list()
        DayEnergyPV <- list()
        BMIEnergyPV <- list()
        DayBMIPV.int <- list()
        DayEnergyPV.int <- list()
        BMIEnergyPV.int <- list()

        #Patients to colors
	colors <- vector()
	cIndex <- 1
	for ( j in 1: nrow(myT))
	{
		if( substr(myT[j,]$Sample.ID,1,1) == "B")
		{
			colors[cIndex] <- "Blue"
		} else if (substr(myT[j,]$Sample.ID,1,1) == "C" )
		{
			colors[cIndex] <- "Red"
		} else
		{
			colors[cIndex] <- "Black"
		}
		cIndex = cIndex + 1
	}

	index <-1

        #Bad hardcode selectors for choosing the patient
        #Patient A
        #myT<- myT[seq(1,67),]
        #Patient B
        #myT<- myT[seq(68,114),]
        #Patient C
        #myT<- myT[seq(115,147),]

        # Should reliably do this even for the mixed case
        myT <- myT[mixedorder(myT[,1]),]

        for(iter in c("A", "B", "C", "ALL"))
           {
        for( i in 2:(ncol(myT) - 14) )
            {
                #Right now this is done too often.
                if(iter == "A")
                {
                    thisT<-myT[seq(1,67),]
                }
                if(iter == "B")
                {
                    thisT<-myT[seq(68,114),]
                }
                if(iter == "C")
                {
                    thisT<-myT[seq(115,147),]
                }
                if(iter == "ALL")
                {
                    thisT <- myT
                }
                #Remove from consideration rare organisms.
		if( sum( thisT[,i] >0 , na.rm=TRUE) > nrow(thisT) /4 )
                    {
                        #Makes these components easier to work with and compute derivative values
                        Day<-thisT$Day
                        ImputedBMI<-thisT$Imputed.BMI
                        EnergyIntake<-thisT$Energy.Intake..kcal.day.
                        taxaType <- as.numeric(thisT[,i])
                        #Compute Shannon diversity and Shannon richness via vegan
                        #myShannon <- diversity(thisT[,i])
                        #myRichness <-
                        if(iter == "ALL")
                        {
                            myFrame <- data.frame(colors, taxaType, Day, ImputedBMI, EnergyIntake)
                        }
                        else
                        {
                            myFrame <- data.frame(taxaType, Day, ImputedBMI, EnergyIntake)
                        }
                        #Should I build an iteration over all of the models too?
                                        #Model building
                        # Simple models without interactions
                        Daylm <- lm(taxaType ~ Day, x = TRUE)
                        BMIlm <- lm(taxaType ~ ImputedBMI, x = TRUE)
                        Energylm <- lm(taxaType ~ EnergyIntake, x = TRUE)

                        DayBMIlm <- lm(taxaType ~ Day + ImputedBMI, x = TRUE)
                        DayEnergylm <- lm(taxaType ~ Day + EnergyIntake, x = TRUE)
                        BMIEnergylm <- lm(taxaType ~ ImputedBMI + EnergyIntake, x = TRUE)

                        myNoInteractLm <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake, x = TRUE)

                                        #Starting interactions
                        DayBMIlm.int <- lm(taxaType ~ Day*ImputedBMI, x = TRUE)
                        DayEnergylm.int <- lm(taxaType ~ Day*EnergyIntake, x = TRUE)
                        BMIEnergylm.int <- lm(taxaType ~ ImputedBMI*EnergyIntake, x = TRUE)

                        #All possible situation
                        fullModel <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake + Day*ImputedBMI + Day*EnergyIntake + ImputedBMI*EnergyIntake + Day*ImputedBMI*EnergyIntake, x = TRUE)
                        #Suggested from fullModel
                        mySuggested <- lm(taxaType ~ Day + ImputedBMI*EnergyIntake, x = TRUE)
                        modelList <- list(Daylm, BMIlm, Energylm, DayBMIlm, DayEnergylm, BMIEnergylm, myNoInteractLm, DayBMIlm.int, DayEnergylm.int, BMIEnergylm.int, fullModel, mySuggested)
                        #With cage/patient effects
                        if(iter == "ALL")
                        {
                            fullModelLME <- lme(taxaType~ Day, method="REML", random = ~1|factor(colors), data = myFrame)
                            modelList <- list(Daylm, BMIlm, Energylm, DayBMIlm, DayEnergylm, BMIEnergylm, myNoInteractLm, DayBMIlm.int, DayEnergylm.int, BMIEnergylm.int, fullModel, mySuggested, fullModelLME)
                        }
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
                        timePValues[index]<-list(summary(Daylm)$coefficients[,4][-1])
                        BMIPValues[index]<-list(summary(BMIlm)$coefficients[,4][-1])
                        EnergyPValues[index]<-list(summary(Energylm)$coefficients[,4][-1])

                        DayBMIPV[index]<-list(summary(DayBMIlm)$coefficients[,4][-1])
                        DayEnergyPV[index]<-list(summary(DayEnergylm)$coefficients[,4][-1])
                        BMIEnergyPV[index]<-list(summary(BMIEnergylm)$coefficients[,4][-1])

                        DayBMIPV.int[index]<-list(summary(DayBMIlm.int)$coefficients[,4][-1])
                        DayEnergyPV.int[index]<-list(summary(DayEnergylm.int)$coefficients[,4][-1])
                        BMIEnergyPV.int[index]<-list(summary(BMIEnergylm.int)$coefficients[,4][-1])

                        noIntPValues[index]<-list(summary(myNoInteractLm)$coefficients[,4][-1])

                        fullPValues[index]<-list(summary(fullModel)$coefficients[,4][-1])

                        names[index] = names(thisT)[i]

                                        #Graphs for each of the models here...
                        for(k in 1:length(modelList))
                            {
                                graphMain = paste(names(thisT)[i],"\n",
                                    paste(names(modelList[[k]]$coefficients[-1]), sep=" "),"\n",
                            "pValues=", paste(format(modelList[[k]]$coefficients[-1], digits=3),sep = " "),
                        sep="")
                                        #Graphing stuff goes here.
                        if (iter == "ALL") {
                            par(mfrow=c(2,1))
                            plot(Day, taxaType, col=colors, pch = 16, main = graphMain)
                            #abline(modelList[[k]])
                            plot(taxaType ~ factor(colors), main=graphMain)
                        }
                        else {
                            plot(Day, taxaType, pch = 16, main=graphMain)
                            #abline(modelList[[k]])
                        }
                    }
                        index = index + 1
		}
            }


        # Building the data.frames to eventually print out the p-values
        DayPV.df <- data.frame(timePValues)
        DayPV.df <- t(DayPV.df)
        dFrame <- data.frame( names, DayPV.df)
#        rownames(dFrame)<- names
	dFrame$adjDay<-  p.adjust( dFrame$Day , method = "BH" )

        BMIPV.df <- data.frame(BMIPValues)
        BMIPV.df <- t(BMIPV.df)
        dFrameBMI <- data.frame(names, BMIPV.df)
 #       rownames(dFrameBMI) <- names
        for (m in 2:dim(dFrameBMI)[2])
        {
            dFrameBMI[,dim(dFrameBMI)[2] + 1] <- p.adjust(dFrameBMI[,m], method = "BH")
            colnames(dFrameBMI)[ncol(dFrameBMI)]<-paste0("adj",colnames(dFrameBMI)[m])
        }

        EnergyPV.df <- data.frame(EnergyPValues)
        EnergyPV.df <- t(EnergyPV.df)
        dFrameEnergy <- data.frame(names, EnergyPV.df)
  #      rownames(dFrameEnergy) <- names
        for (m in 2:dim(dFrameEnergy)[2])
        {
            dFrameEnergy[,dim(dFrameEnergy)[2] + 1] <- p.adjust(dFrameEnergy[,m], method = "BH")
            colnames(dFrameEnergy)[ncol(dFrameEnergy)]<-paste0("adj",colnames(dFrameEnergy)[m])
        }

        noIntPV.df <- data.frame(noIntPValues)
        noIntPV.df <- t(noIntPV.df)
        dFrameNoInt <- data.frame(names, noIntPV.df)
   #     rownames(dFrameNoInt)<- names
        for (m in 2:dim(dFrameNoInt)[2])
        {
            dFrameNoInt[,dim(dFrameNoInt)[2] + 1] <- p.adjust(dFrameNoInt[,m], method = "BH")
            colnames(dFrameNoInt)[ncol(dFrameNoInt)]<-paste0("adj",colnames(dFrameNoInt)[m])
        }

        fullPV.df <- data.frame(fullPValues)
        fullPV.df <- t(fullPV.df)
        dFrameFull <- data.frame(names, fullPV.df)
    #    rownames(dFrameFull) <- names
        for (m in 2:dim(dFrameFull)[2])
        {
           dFrameFull[,dim(dFrameFull)[2] + 1] <- p.adjust(dFrameFull[,m], method = "BH")
           colnames(dFrameFull)[ncol(dFrameFull)]<-paste0("adj",colnames(dFrameFull)[m])
        }

                                        #Finally, writing out the p-values and BH adjusted p-values
        #Still need to add in the other cases
	write.table(dFrame, file= paste( "pValuesLongitudinalDayModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameBMI, file= paste( "pValuesLongBMIModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameEnergy, file= paste( "pValuesLongEnergyModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        #write.table(BMI, file= paste( "pValuesLongitudinalDayModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        #write.table(dFrameBMI, file= paste( "pValuesLongBMIModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        #write.table(dFrameEnergy, file= paste( "pValuesLongEnergyModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameNoInt, file = paste( "pValuesLongNoInt_", iter, "_", t, ".txt",sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameFull, file = paste("pValuesLongFullModel_", iter, "_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

        #TODO Do cutoffs at the 0.05 level and print out those reduced files as well.

    }
    dev.off()
}
