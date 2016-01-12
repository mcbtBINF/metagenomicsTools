rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")

library("Kendall")
library("vegan")
library("nlme")

taxaLevels <- c("phylum","class","order","family","genus")
#Right Now I'm being lazy with the changing of the variable names

for(t in taxaLevels )
{
        t<-taxaLevels[[1]]
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
        myT<- myT[seq(1,67),]

                                        #Patient B
#        myT<- myT[seq(68,114),]
                                        #Patient C
                                        #        myT<- myT[seq(115,147),]
#Excludes other measurement columns
        pval.list <- list()
        rsquared.list <- list()

        for( i in 2:(ncol(myT) - 14) )
            {
                                        #Remove from consideration rare organisms.

		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                                                Day<-myT$Day
                        ImputedBMI<-myT$Imputed.BMI
                                                EnergyIntake<-myT$Energy.Intake..kcal.day.
                                                taxaType <- myT[,i]

                                                # Simple model which doesn't explain much of the variance
                                                myLm <- lm(taxaType ~ Day, x = TRUE)


                                                pval.list[[i]]<-summary(myLm)$coefficients[,4]
                                                rsquared.list[[i]]<-summary(myLm)$r.squared

                                                anova(myLm)
                                                plot(Day, taxaType)
                                                abline(myLm)
                                        #Simple No Interaction
                                                BMILm <- lm(taxaType ~ ImputedBMI + EnergyIntake, x = TRUE)
                                                myNoInteractLm <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake, x = TRUE)

                                                #All possible situation
                                                fullModel <- lm(taxaType ~ Day + ImputedBMI + EnergyIntake + Day*ImputedBMI + Day*EnergyIntake + ImputedBMI*EnergyIntake, x=TRUE)
                                        #Suggested from fullModel
                                                mySuggested <- lm(taxaType ~ Day + ImputedBMI*EnergyIntake, x = TRUE)
                                        #Colors as a proxy for patient.
                                        # Working on single color at this point
#                        print(colors[[i]])


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
#	dFrame$adjPatient<-  p.adjust( dFrame$patientPValues, method = "BH" )
#      	dFrame <- data.frame( names,patientPValues, BMIPValues, interactionPValues )
 #      dFrame <- dFrame [order(dFrame$BMIPValues),]
#	dFrame$adjBMI<-  p.adjust( dFrame$BMIPValues , method = "BH" )
#	dFrame$adjPatient<-  p.adjust( dFrame$patientPValues, method = "BH" )


	write.table( file= paste( "pValuesLongitudinalModel_", t, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")

#Temporarily just looking at one taxonomic level
}
