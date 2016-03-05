rm(list=ls())

# setwd("C:\\ianLongitudinal")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrolData/Carroll_Longitudinal")
library("Kendall")

taxaLevels <- c("phylum","class","order","family","genus")


#t <- "family"
for(t in taxaLevels )
{
	pdf(paste(t,"_plots.pdf", sep="") )
#	inFileName <- paste(t,"LogNormalwithMetadata.txt", sep="")
      	inFileName <- paste(t,"LogNormalwithMetadata_Edit.txt", sep="")
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

	for( i in 3:ncol(myT))
	{
		if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
		{
			myLm <- lm( myT[,i] ~ colors *  myT$Day )

			myAnova <- anova(myLm)

			patientPValues[index] <- myAnova$"Pr(>F)"[1]
			timePValues[index] <- myAnova$"Pr(>F)"[2]
			interactionPValues[index] <- myAnova$"Pr(>F)"[3]
			names[index] <- names(myT)[i]

			myLabel <- paste(names(myT)[i] , "\n", "p patient = " ,
			format(patientPValues[index],digits=3) ,
					"\n p time = " , format(timePValues[index],digits=3), "\n", "p interaction = " ,
					format(interactionPValues[index],digits=3))

			plot(myT[,i] ~ myT$Day, main = myLabel, col=colors)

                        coefs <- coef(myLm)
                        abline( a=coefs[1] + coefs[2], b=coefs[5], col="BLUE")
                        abline( a=coefs[1] + coefs[3], b=coefs[6], col="RED")
                        index = index + 1

		}
	}
	dev.off()

	dFrame <- data.frame( names,patientPValues, timePValues ,interactionPValues )
        dFrame <- dFrame [order(dFrame$timePValues),]
	dFrame$adjTime<-  p.adjust( dFrame$timePValues , method = "BH" )
	dFrame$adjPatent<-  p.adjust( dFrame$patientPValues, method = "BH" )

	write.table( file= paste( "pValuesLongitudinalModel_", t, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}
