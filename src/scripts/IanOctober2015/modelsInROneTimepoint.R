rm(list=ls())
library("lmtest")
library("nlme")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/BombCalorimetry/")

taxaLevels <- c("phylum","class","order","family","genus")

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa,  "_paired_metadata.txt", sep ="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",1), rep("numeric", numCols-1))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

	myT <- myT[ myT$Time == "1" &  ! is.na(myT$calorimetryData), ]

	names <- vector()
	pValuesSubject <- vector()
	pValuesCalorimetry <- vector()
	meanBug <- vector()
	index <- 1
	pdf( paste(taxa, "plotsTimepoint1.pdf", sep=""))

	for( i in 6:numCols)
		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- log10( myT[,i] + 0.00001)
			meanBug[index] <- mean(bug)
			time <- factor(myT$timepoint)
			patientID <- myT$patientID
			calorimetry<- myT$calorimetryData

			myFrame <- data.frame(bug, time, patientID, calorimetry)

			fullModel <- lm( bug~  calorimetry)

			pValuesCalromitery[index] <- anova(fullModel)$"Pr(>F)"[1]
			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i],
				" pValuesCalorimetry= ", format(pValuesCalromitery[index],digits=3))

			plot( bug ~ calorimetry, ylab = names[index],
					main = graphMain )
			index=index+1

		}

	dFrame <- data.frame( names, pValuesCalromitery ,meanBug)
	dFrame <- dFrame [order(dFrame$pValuesCalromitery),]
	dFrame$adjustedpValuesCalromitery <- p.adjust( dFrame$pValuesCalromitery, method = "BH" )
	write.table(dFrame, file=paste("pValuesFor", taxa, "OnlyTime1.txt",sep=""), sep="\t",row.names=FALSE)
		dev.off()
}
