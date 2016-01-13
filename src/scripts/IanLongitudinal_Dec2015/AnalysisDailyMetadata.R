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
#      t<-taxaLevels[[1]]
#	pdf(paste(t,"_plots.pdf", sep="") )
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

        #Patients to colors
	colors <- vector()
	cIndex <-1
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
                                        #Bad hardcode selectors
        #Choosing a patient
                                        #Patient A
#        myT<- myT[seq(1,67),]

                                        #Patient B
#        myT<- myT[seq(68,114),]
                                        #Patient C
#        myT<- myT[seq(115,147),]
                                        #Excludes other measurement columns
        # Should reliably do this even for the mixed case
        myT <- myT[mixedorder(myT[,1]),]

        for( i in 2:(ncol(myT) - 14) )
            {
                                        #Remove from consideration rare organisms.

		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                                                Day<-myT$Day
                        ImputedBMI<-myT$Imputed.BMI
                                                EnergyIntake<-myT$Energy.Intake..kcal.day.
                                                taxaType <- as.numeric(myT[,i])

                                                myFrame <- data.frame(colors, taxaType, Day, ImputedBMI, EnergyIntake)
                                             #   if(modelNum == 0) {
                                        # Simple model which doesn't explain much of the variance

                                                fullModelLME <- lme(taxaType~ Day, method="REML", random = ~1|factor(colors), data = myFrame)

                                                myLm <- lm(taxaType ~ Day, x=TRUE)
                                            #}
       #                                                pval.list[[i]]<-summary(myLm)$coefficients[,4]
 #                                               rsquared.list[[i]]<-summary(myLm)$r.squared

                                             #   anova(myLm)
                                              #  plot(Day, taxaType)
                                               # abline(myLm)
                                        #Simple No Interaction
                                                BMIlm <- lm(taxaType ~ ImputedBMI, x = TRUE)
                                                Energylm <- lm(taxaType ~ EnergyIntake, x=TRUE)
                                                myNoInteractLm <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake, x = TRUE)

                                        #All possible situation
                                               # if(modelNum == 21) {
                                                    fullModel <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake + Day*ImputedBMI + Day*EnergyIntake + ImputedBMI*EnergyIntake + Day*ImputedBMI*EnergyIntake, x = TRUE)
                                               # }
                                        #Suggested from fullModel
                                                mySuggested <- lm(taxaType ~ Day + ImputedBMI*EnergyIntake, x = TRUE)

                                        # Second corner
     #                                           time2ndOrderPseud <- lm( myT$LogPseudo ~ poly( myT$sampleDays,2), x=TRUE)
                                                # This is probably just developed for the simple regular intervals of time points.

                                        # Autocorrelation Code
      #                                          M1 <- gls( taxaType ~ EnergyIntake , data= myFrame )
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


                                                timePValues[index]<-summary(myLm)$coefficients[,4][[2]]
                                                BMIPValues[index]<-list(summary(BMIlm)$coefficients[,4][-1])
EnergyPValues[index]<-list(summary(Energylm)$coefficients[,4][-1])

noIntPValues[index]<-list(summary(myNoInteractLm)$coefficients[,4][-1])
                                                fullPValues[index]<-list(summary(fullModel)$coefficients[,4][-1])

                                        #print(summary(myLm)$coefficients[,4][[2]])

                                        #                                                timePValues[index] <- anova(myLm)$"p-value"[2]

                                                names[index] = names(myT)[i]


                                                # Graphing stuff goes here.

                                                index = index + 1

                                        #Compute Shannosn diversity and Shannon richness via vegan
#                        myShannon <- diversity(myT[,i])
#                        myRichness <-

                                        #myLm <- lm( myT[,i] ~ colors *  myT$Day )
#                        myLm <- lm( myT[,i] ~  myT$Day )
#                        checkAll <- lm(myT[,i] ~ #colors +
 #                                          myT$Day + myT$Imputed.BMI + myT$Energy.Intake..kcal.day. +
                                        # Interactions
                                           #colors*myT$Day + colors*myT$Imputed.BMI + colors*myT$Energy.Intake..kcal.day. +
  #                                             myT$Day*myT$Imputed.BMI + myT$Imputed.BMI*myT$Energy.Intake..kcal.day. + myT$Day*myT$Energy.Intake..kcal.day. +
   #                                                myT$Day*myT$Imputed.BMI*myT$Energy.Intake..kcal.day.
    #                                   )
     #                                           rescheckAll <- rstandard(checkAll)
                                        # Plot the fitted versus the standardized residuals
#                                                pdf( paste( paste(t, names(myT)[i], sep="_"),"resid.pdf", sep="_"))
 #                                               par(mfrow=c(2,2))

  #                                                  plot(fitted(checkAll), rescheckAll)
                                        # Plot each explanatory variable versus the standardized residuals

   #                                             plot(Day, rescheckAll)
    #                                            plot(ImputedBMI, rescheckAll)
     #                                           plot(EnergyIntake, rescheckAll)

      #                                          dev.off()

#                        checkAllgls <- gls(taxaType ~ #colors +
#                                           Day + ImputedBMI + EnergyIntake +
                                        # Interactions
                                           #colors*myT$Day + colors*myT$Imputed.BMI + colors*myT$Energy.Intake..kcal.day. +
#                                               Day*ImputedBMI + ImputedBMI*EnergyIntake + Day*EnergyIntake +
#                                                  Day*ImputedBMI*EnergyIntake
#                                      )


#			myAnova <- anova(myLm)

#			patientPValues[index] <- myAnova$"Pr(>F)"[1]
 #                   timePValues[index] <- myAnova$"Pr(>F)"[2]
                    #BMIPValues[index] <- myAnova$"Pr(>F)"[2]
#			interactionPValues[index] <- myAnova$"Pr(>F)"[3]
#			names[index] <- names(myT)[i]

#			myLabel <- paste(names(myT)[i] , "\n", "p patient = " ,
#			format(patientPValues[index],digits=3) ,
#                                                                                 "\n p time = " , format(timePValues[index],digits=3), "\n",
                   #                      "\n p BMI = ",
                    #                     format(BMIPValues[index],digits=3), "\n",
 #                                        "p interaction = " ,
#					format(interactionPValues[index],digits=3))

 #                  plot(myT[,i] ~ myT$Day, main = myLabel, col=colors, ylim=c(0,5))
                                        #plot(myT[,i] ~ myT$BMI, main = myLabel, col=colors)
                   #plot(myT$BMI ~ myT$Day, main = myLabel, col=colors)


 #                       coefs <- coef(myLm)
                     #   abline( a=coefs[1] + coefs[2], b=coefs[5], col="BLUE")
                     #   abline( a=coefs[1] + coefs[3], b=coefs[6], col="RED")
 #                       index = index + 1

		}
            }
 #                   }

#	dev.off()


	dFrame <- data.frame( names, timePValues) # ,interactionPValues )
        dFrame <- dFrame [order(dFrame$timePValues),]
	dFrame$adjTime<-  p.adjust( dFrame$timePValues , method = "BH" )

        BMIPV.df <- data.frame(BMIPValues)
        BMIPV.df <- t(BMIPV.df)
        dFrameBMI <- data.frame(names, BMIPV.df)
        rownames(dFrameBMI) <- names
        for (m in 2:dim(dFrameBMI)[2]) {
            dFrameBMI[,dim(dFrameBMI)[2] + 1] <- p.adjust(dFrameBMI[,m], method = "BH")
            colnames(dFrameBMI)[ncol(dFrameBMI)]<-paste0("adj",colnames(dFrameBMI)[m])
        }

        EnergyPV.df <- data.frame(EnergyPValues)
        EnergyPV.df <- t(EnergyPV.df)
        dFrameEnergy <- data.frame(names, EnergyPV.df)
        rownames(dFrameEnergy) <- names
        for (m in 2:dim(dFrameEnergy)[2]) {
            dFrameEnergy[,dim(dFrameEnergy)[2] + 1] <- p.adjust(dFrameEnergy[,m], method = "BH")
            colnames(dFrameEnergy)[ncol(dFrameEnergy)]<-paste0("adj",colnames(dFrameEnergy)[m])
        }


        noIntPV.df <- data.frame(noIntPValues)
        noIntPV.df <- t(noIntPV.df)
        dFrameNoInt <- data.frame(names, noIntPV.df)
        rownames(dFrameNoInt)<- names
        for (m in 2:dim(dFrameNoInt)[2]) {
            dFrameNoInt[,dim(dFrameNoInt)[2] + 1] <- p.adjust(dFrameNoInt[,m], method = "BH")
            colnames(dFrameNoInt)[ncol(dFrameNoInt)]<-paste0("adj",colnames(dFrameNoInt)[m])
        }

        fullPV.df <- data.frame(fullPValues)
        fullPV.df <- t(fullPV.df)
        dFrameFull <- data.frame(names, fullPV.df)
        rownames(dFrameFull) <- names
        for (m in 2:dim(dFrameFull)[2]) {
            dFrameFull[,dim(dFrameFull)[2] + 1] <- p.adjust(dFrameFull[,m], method = "BH")
           colnames(dFrameFull)[ncol(dFrameFull)]<-paste0("adj",colnames(dFrameFull)[m])
        }



                                        #	dFrame$adjPatient<-  p.adjust( dFrame$patientPValues, method = "BH" )
#      	dFrame <- data.frame( names,patientPValues, BMIPValues, interactionPValues )
 #      dFrame <- dFrame [order(dFrame$BMIPValues),]
#	dFrame$adjBMI<-  p.adjust( dFrame$BMIPValues , method = "BH" )
#	dFrame$adjPatient<-  p.adjust( dFrame$patientPValues, method = "BH" )


	write.table(dFrame, file= paste( "pValuesLongitudinalSimpleModel_ALL_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameBMI, file= paste( "pValuesLongBMIModel_ALL_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameEnergy, file= paste( "pValuesLongEnergyModel_ALL_", t, ".txt", sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameNoInt, file = paste( "pValuesLongNoInt_ALL_",t, ".txt",sep=""), row.names=FALSE, sep="\t")
        write.table(dFrameFull, file = paste("pValuesLongFullModel_ALL_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

#Temporarily just looking at one taxonomic level
}
