## Need to add the largest technical replicate selection
## Need to add the nice plotting

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/qiime/closed_reference_merged/countByTaxonomy")


getColumnIndex <- function(myT, s)
{
    myNames <- names(myT)

    for( i in 1: length(myNames))
        if ( myNames[i] == s)
            return (i)
}


## i <- 2
for( i in 2:7 )
{
    fileName <- paste("LogNormwithMetadata_L_", i, ".txt", sep = "")
    myT <- read.table(fileName, sep="\t", header=TRUE)

    ## This may cause problems with the data-types of the columns
    myT <- data.frame(lapply(myT, function(x) {
        gsub("ND|<LOQ", 0, x)
    }))

    ## pValueLocationsFromFull <- vector()
    ## pValueUpDownFromFull<- vector()
    ## pValuesTimepointFromFull <- vector()
    ## pValuesAntibioticFromFull <- vector()
    ## pValuesUpDownLocationInteraction <- vector()

    metadataStart <- getColumnIndex(myT,"sequenceCount")

    ## pValuesLocation <- vector()
    ## names <- vector()
    ## index <-1
    ## pdf( paste("Antibiotic_AssociationL", i, ".pdf"))

    ## Get the technical replicate from the deepest sequenced sample.
    savemyT <- myT[FALSE,]
    for (eachVal in unique(myT$Sample.ID)) {
        myTperSample <- myT[myT$Sample.ID == eachVal,]
        rep<-myTperSample[which(as.numeric(levels(myTperSample$sequenceCount))[myTperSample$sequenceCount] == max(as.numeric(levels(myTperSample$sequenceCount)[myTperSample$sequenceCount]))),]
        savemyT <- rbind(savemyT, rep)
    }
    myT <- savemyT

    justFirst <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 1))
    justFirst<-gsub("/16", "/2016", justFirst)
    justFirst <- as.numeric(as.POSIXlt(strptime(justFirst, "%m/%d/%Y")))
    myT$Timestamp <- justFirst

    ## Go over the 10 antibiotics
    for( k in (metadataStart + 16):(metadataStart + 25))
    {
        pValuesLocation <- vector()
        names <- vector()

        pValueLocationsFromFull <- vector()
        pValueUpDownFromFull<- vector()
        pValuesTimepointFromFull <- vector()
        pValuesAntibioticFromFull <- vector()
        pValuesUpDownLocationInteraction <- vector()

        index <-1
        pdf( paste(names(myT)[k], "Antibiotic_Association_L", i, ".pdf"))

	for( j in 2: (metadataStart -1) )
	{
            par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
            bug <- as.numeric(levels(myT[,j]))[myT[,j]]
            if( sum(bug != 0) > nrow(myT) / 4)
            {
                justStreams <- myT[(myT$Location == "Mallard Creek" |myT$Location ==  "Sugar Creek") &
                                   (myT$Sample == "UP A" |myT$Sample == "UP B" |
                                    myT$Sample == "DS A" | myT$Sample == "DS B"), ]

                locations <- factor(justStreams$Location)
                streamBugs <- as.numeric(levels(justStreams[,j]))[justStreams[,j]]
                thisAntibiotic <- as.numeric(levels(justStreams[,k]))[justStreams[,k]]
                ## thisAntibiotic <- as.numeric(justStreams[,k])
                ## myLm <- lm( streamBugs ~ locations)
                names[index] <- names(justStreams)[j]
                pValuesLocation[index] <- wilcox.test( streamBugs[myT$Location == "Mallard Creek"] ,
                                                      streamBugs[myT$Location == "Sugar Creek"])$p.value
                ## mainText = paste(names[index], format(pValuesLocation[index], digits=3))
                boxplot( streamBugs ~ locations, main = paste("Stream (non-parametric) \nuncorrected p-value",  format(pValuesLocation[index], digits=3)))

                updownBinary <- ifelse( justStreams$Sample == "DS A" | justStreams$Sample == "DS B" ,
                                       "down", "up"  )

                myFrame <- data.frame(streamBugs,locations, updownBinary, justStreams$Timestamp)

                stripchart(streamBugs ~ locations,
                           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

                fullModel <- lm( streamBugs ~ locations * updownBinary + justStreams$Timestamp + thisAntibiotic)
                anAnova <- anova(fullModel)

                pValueLocationsFromFull[index] <- anAnova$"Pr(>F)"[1]
                pValueUpDownFromFull[index] <- anAnova$"Pr(>F)"[2]
                pValuesTimepointFromFull[index] <- anAnova$"Pr(>F)"[3]
                pValuesAntibioticFromFull[index] <- anAnova$"Pr(>F)"[4]
                pValuesUpDownLocationInteraction[index] <- anAnova$"Pr(>F)"[5]

                boxplot( streamBugs ~ updownBinary, main = paste("Up/Down \nuncorrected p-value",  format(pValueUpDownFromFull[index], digits=3)))
                stripchart(streamBugs ~ updownBinary,
                           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
                ## boxplot( streamBugs ~ as.POSIXct(justStreams$Timestamp, origin = "1970-01-01"), main = format(pValuesTimepointFromFull[index], digits=3))
                ## stripchart(streamBugs ~ as.POSIXct(justStreams$Timestamp, origin = "1970-01-01"),
                ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
                boxplot( streamBugs ~ justStreams$Timepoint, main = paste("Time \nuncorrected p-value", format(pValuesTimepointFromFull[index], digits=3)))
                stripchart(streamBugs ~ justStreams$Timepoint,
                           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
                ## plot(justStreams$Timepoint, streamBugs, col=ifelse(justStreams$Location == "Mallard Creek", "red", "blue"), pch = ifelse(updownBinary == "up", 24, 25), main="Red=Mallard \npointed-up triangle=upstream")
                boxplot( streamBugs ~ thisAntibiotic, main = paste("Antibiotic ", names(justStreams)[k], "\nuncorrected p-value", format(pValuesTimepointFromFull[index], digits=3)))
                stripchart(streamBugs ~ thisAntibiotic,
                           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])


                ## mtext(unlist(strsplit(names[index],split="\\."))[i], outer=TRUE, cex = 1.5)
                mtext(names[index], outer=TRUE, cex = 0.75)

                index = index + 1
            }
	}

	hist(pValuesLocation,breaks=20, main="Non-parametric Stream")

	hist(pValueLocationsFromFull,breaks=20, main="Stream")
	hist(pValueUpDownFromFull,breaks=20, main="Up/Downstream")
	hist(pValuesTimepointFromFull,breaks=20, main="Timepoint")
      	hist(pValuesAntibioticFromFull,breaks=20, main=paste("Antibiotic Concentration", names(justStreams)[k]))
	hist(pValuesUpDownLocationInteraction,breaks=20, main="Stream and Up/Down")

	dev.off()
	myFrame <- data.frame(names, pValuesLocation,                              pValueLocationsFromFull,pValueUpDownFromFull,pValuesTimepointFromFull, pValuesAntibioticFromFull,
                              pValuesUpDownLocationInteraction)
	myFrame <- myFrame [order(myFrame$pValuesLocation),]
	myFrame$adjustedPValuesLocation <- p.adjust( myFrame$pValuesLocation, method = "BH" )

	myFrame$adjustedpValueLocationsFromFull <- p.adjust( myFrame$pValueLocationsFromFull, method = "BH" )
	myFrame$adjustedpValueUpDownFromFull <- p.adjust( myFrame$pValueUpDownFromFull, method = "BH" )
	myFrame$adjustedpValuesTimepointFromFull <- p.adjust( myFrame$pValuesTimepointFromFull, method = "BH" )
       	myFrame$adjustedpValuesAntibioticFromFull <- p.adjust( myFrame$pValuesAntibioticFromFull, method = "BH" )
	myFrame$adjustedpValuesUpDownLocationInteraction <- p.adjust( myFrame$pValuesUpDownLocationInteraction, method = "BH" )

	write.table(myFrame, file=paste("L_", i, "_", names(myT)[k], "_AntibioticAssociation_pValues.txt",sep=""), sep="\t",row.names=FALSE)
    }
}

