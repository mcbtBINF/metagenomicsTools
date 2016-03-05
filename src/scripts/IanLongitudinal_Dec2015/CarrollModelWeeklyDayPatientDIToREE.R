rm(list=ls())

# setwd("C:\\ianLongitudinal")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")
library("Kendall")

taxaLevels <- c("phylum","class","order","family","genus")


#t <- "family"
for(t in taxaLevels )
{
	pdf(paste(t,"_weekly_DayPatientActiveEE_plots.pdf", sep="") )
#	inFileName <- paste(t,"LogNormalwithMetadata.txt", sep="")
      	inFileName <- paste(t,"LogNormalwithMetadataWeekly_NearestSamplewithDay_Edit.txt", sep="")
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

	for( i in 3:(ncol(myT) - 14))
	{
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
                    {
                        DIT <- myT$DIT..kcal.day.
                        REE <- myT$REE..kcal.day.
                        ActiveEE <- myT$Active.EE..kcal.
                        DITREE <- DIT / REE
                        myLm <- lm( myT[,i] ~ myT$Day*patient*ActiveEE)

                        myLmpVal[index] <- list(summary(myLm)$coefficients[,4][-1])
			myAnova <- anova(myLm)

#			patientPValues[index] <- myAnova$"Pr(>F)"[1]
#			timePValues[index] <- myAnova$"Pr(>F)"[2]
 #                   REEPValues[index] <- myAnova$"Pr(>F)"[3]
                    # 4 if only first two interacting, 7 if all three interact
  #                      interactionPValues[index] <- myAnova$"Pr(>F)"[4]
			names[index] <- names(myT)[i]

#			myLabel <- paste(names(myT)[i] , "\n", "p patient = " ,
			#format(patientPValues[index],digits=3) ,
#                                         " p time = " , format(timePValues[index],digits=3),
#                                         "\n p ActiveEE_kcal = ", format(REEPValues[index], digits=3),
#                                         "\n", "p interaction = " ,
#					format(interactionPValues[index],digits=3))


#                     plot(myT$Day, myT[,i], main = graphMain, col=colors)

 #                    coefs <- coef(myLm)
                     index = index + 1

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

	write.table(dFramemyLm, file= paste( "pValuesLongitudinalModelWeekly_DayPatientActiveEE_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

#        dev.off()
}
