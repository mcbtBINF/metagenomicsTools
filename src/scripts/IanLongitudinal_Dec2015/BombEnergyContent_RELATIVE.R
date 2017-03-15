rm(list=ls())

## setwd("C:\\ianLongitudinal")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")
library("Kendall")
library("vegan")
library("lmtest")
library("pscl")
library("nlme")
library("gtools")

taxaLevels <- c("phylum","class","order","family","genus")

## t <- "family"
for(t in taxaLevels )
{

    ## inFileName <- paste(t,"LogNormalwithMetadata.txt", sep="")
    inFileName <- paste(t,"_LogNormalwithMetadataWeekly_NearestSamplewithDay_withCals.txt", sep="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    myColClasses <- c("character", rep("numeric", numCols-1))
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
    ##        myT <- myT[is.na(myT$Day) == FALSE,]
    ##        myT <- myT[mixedorder(myT[,1]),]
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

    pdf(paste(t,"_weekly_PatientEnergyContent_RELATIVE_plots.pdf", sep="") )

    for( i in 2:(ncol(myT) - 15))
        {
            ## I know that this should change, but I'm not sure how it should change.
            if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                {
                    DIT <- myT$DIT..kcal.day.
                    REE <- myT$REE..kcal.day.
                    ActiveEE <- myT$Active.EE..kcal.
                    DITREE <- DIT / REE
                    EnergyContent <- myT$Energy.Content..cal.g/myT$Energy.Intake..kcal.day.
                    ## EnergyContent <- myT$
                    myLm <- lm( myT[,i] ~ patient*EnergyContent)

                    myLmpVal[index] <- list(summary(myLm)$coefficients[,4][-1])
                    myAnova <- anova(myLm)

                    ## patientPValues[index] <- myAnova$"Pr(>F)"[1]
                    ## timePValues[index] <- myAnova$"Pr(>F)"[2]
                    ## REEPValues[index] <- myAnova$"Pr(>F)"[3]
                    ## 4 if only first two interacting, 7 if all three interact
                    ## interactionPValues[index] <- myAnova$"Pr(>F)"[4]
                    names[index] <- names(myT)[i]

                    ## myLabel <- paste(names(myT)[i] , "\n", "p patient = " ,
                    ## format(patientPValues[index],digits=3) ,
                    ## " p time = " , format(timePValues[index],digits=3),
                    ## "\n p ActiveEE_kcal = ", format(REEPValues[index], digits=3),
                    ## "\n", "p interaction = " ,
                    ## format(interactionPValues[index],digits=3))

                    meanA <- mean(as.numeric(myT[colors == "Black", i]))
                    meanB <- mean(as.numeric(myT[colors == "Blue", i]))
                    meanC <- mean(as.numeric(myT[colors == "Red", i]))

                    graphMain <- names(myT)[i]
                    ## graphMain <- paste0(names(myT)[i],"\n", "meanA=", meanA, "\n", "meanB=", meanB, "\n", "meanC=", meanC)#, "\n", "coef", myLm$coef)
                    plot(EnergyContent, myT[,i], main = graphMain, col=colors, xlab="Relative Energy Content", ylab = "Log-Abundance", pch=16, ylim=c(0,5))

                    ## abline(a = meanA
                    ##        ## myLm$coef[1]
                    ##      , b = myLm$coef[4])
                    ## abline(a = meanB
                    ##        ## myLm$coef[1] + myLm$coef[2] #+ myLm$coef[6]
                    ##      , b = myLm$coef[4] + myLm$coef[5] #+ myLm$coef[2]
                    ##      , col="BLUE")
                    ## abline(a = meanC
                    ##        ## myLm$coef[1] + myLm$coef[3] #+ myLm$coef[4]
                    ##      , b = myLm$coef[4] + myLm$coef[6] #+ myLm$coef[2]
                    ##      , col="RED")

                    ## coefs <- coef(myLm)
                    index = index + 1
                    ## dev.off()
		}
	}
    dev.off()

    dFrame <- data.frame( myLmpVal)
    dFrame <- t(dFrame)
    dFramemyLm <- data.frame(names, dFrame)

    for (m in 2:dim(dFramemyLm)[2])
        {
            dFramemyLm[,dim(dFramemyLm)[2] + 1] <- p.adjust(dFramemyLm[,m], method = "BH")
            colnames(dFramemyLm)[ncol(dFramemyLm)]<-paste0("adj", colnames(dFramemyLm)[m])
        }

    write.table(dFramemyLm, file= paste( "pValuesLongitudinalModelWeekly_PatientEnergyContent_RELATIVE_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

    ## dev.off()
}
