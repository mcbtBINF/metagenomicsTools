## Need to add the largest technical replicate selection
## Need to add the nice plotting

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/qiime/de_novo_merged/denovoLevels")


getColumnIndex <- function(myT, s)
{
	myNames <- names(myT)

	for( i in 1: length(myNames))
		if ( myNames[i] == s)
			return (i)
}


i <- 2
	fileName <- paste("denovoLogNormwithMetadata_L_", i, ".txt", sep = "")
	myT <- read.table(fileName, sep="\t", header=TRUE)

msCorrect <- read.table("LODLOQcorrection.txt", sep="\t", header=TRUE)
## This will be changed
## The data is on the drive in file: UNCC161_LOD.xlsx

myT$Ertapenem..ng.L. <- gsub("ND", msCorrect$Ertapenem[1], myT$Ertapenem..ng.L.)
myT$Ertapenem..ng.L. <- gsub("<LOQ", msCorrect$Ertapenem[2], myT$Ertapenem..ng.L.)
myT$Ertapenem..ng.L. <- as.numeric(myT$Ertapenem..ng.L.)
myT$Amoxicillin..ng.L. <- gsub("ND", msCorrect$Amoxicillin[1], myT$Amoxicillin..ng.L.)
myT$Amoxicillin..ng.L. <- gsub("<LOQ", msCorrect$Amoxicillin[2], myT$Amoxicillin..ng.L.)
myT$Amoxicillin..ng.L. <- as.numeric(myT$Amoxicillin..ng.L.)
myT$Ciprofloxacin..ng.L. <- gsub("ND", msCorrect$Ciprofloxacin[1], myT$Ciprofloxacin..ng.L.)
myT$Ciprofloxacin..ng.L. <- gsub("<LOQ", msCorrect$Ciprofloxacin[2], myT$Ciprofloxacin..ng.L.)
myT$Ciprofloxacin..ng.L. <- as.numeric(myT$Ciprofloxacin..ng.L.)
myT$Doxycycline..ng.L. <- gsub("ND", msCorrect$Doxycycline[1], myT$Doxycycline..ng.L.)
myT$Doxycycline..ng.L. <- gsub("<LOQ", msCorrect$Doxycycline[2], myT$Doxycycline..ng.L.)
myT$Doxycycline..ng.L. <- as.numeric(myT$Doxycycline..ng.L.)
myT$Azithromycin..ng.L. <- gsub("ND", msCorrect$Azithromycin[1], myT$Azithromycin..ng.L.)
myT$Azithromycin..ng.L. <- gsub("<LOQ", msCorrect$Azithromycin[2], myT$Azithromycin..ng.L.)
myT$Azithromycin..ng.L. <- as.numeric(myT$Azithromycin..ng.L.)
myT$Clindamycin..ng.L. <- gsub("ND", msCorrect$Clindamycin[1], myT$Clindamycin..ng.L.)
myT$Clindamycin..ng.L. <- gsub("<LOQ", msCorrect$Clindamycin[2], myT$Clindamycin..ng.L.)
myT$Clindamycin..ng.L. <- as.numeric(myT$Clindamycin..ng.L.)
myT$Sulfamethoxazole..ng.L. <- gsub("ND", msCorrect$Sulfamethoxazole[1], myT$Sulfamethoxazole..ng.L.)
myT$Sulfamethoxazole..ng.L. <- gsub("<LOQ", msCorrect$Sulfamethoxazole[2], myT$Sulfamethoxazole..ng.L.)
myT$Sulfamethoxazole..ng.L. <- as.numeric(myT$Sulfamethoxazole..ng.L.)
myT$Cephalexin..ng.L. <- gsub("ND", msCorrect$Cephalexin[1], myT$Cephalexin..ng.L.)
myT$Cephalexin..ng.L. <- gsub("<LOQ", msCorrect$Cephalexin[2], myT$Cephalexin..ng.L.)
myT$Cephalexin..ng.L. <- as.numeric(myT$Cephalexin..ng.L.)
myT$Trimethoprim..ng.L. <- gsub("ND", msCorrect$Trimethoprim[1], myT$Trimethoprim..ng.L.)
myT$Trimethoprim..ng.L. <- gsub("<LOQ", msCorrect$Trimethoprim[2], myT$Trimethoprim..ng.L.)
myT$Trimethoprim..ng.L. <- as.numeric(myT$Trimethoprim..ng.L.)
myT$Levofloxacin..ng.L. <- gsub("ND", msCorrect$Levofloxacin[1], myT$Levofloxacin..ng.L.)
myT$Levofloxacin..ng.L. <- gsub("<LOQ", msCorrect$Levofloxacin[2], myT$Levofloxacin..ng.L.)
myT$Levofloxacin..ng.L. <- as.numeric(myT$Levofloxacin..ng.L.)

## myT <- data.frame(lapply(myT, function(x) {
##             gsub("ND|<LOQ", 0, x)
##         }))

	pValueLocationsFromFull <- vector()
	pValueUpDownFromFull<- vector()
	pValuesTimepointFromFull <- vector()
	pValuesUpDownLocationInteraction <- vector()

	metadataStart <- getColumnIndex(myT,"sequenceCount")

	pValuesLocation <- vector()
	names <- vector()
	index <-1
	pdf( paste("denovoAntibiotics_L", i, ".pdf"))

        ## Get the technical replicate from the deepest sequenced sample.
        savemyT <- myT[FALSE,]
        for (eachVal in unique(myT$Sample.ID)) {
            myTperSample <- myT[myT$Sample.ID == eachVal,]
            rep<-myTperSample[which(as.numeric(myTperSample$sequenceCount) == max(as.numeric(myTperSample$sequenceCount))),]
            savemyT <- rbind(savemyT, rep)
        }
        myT <- savemyT

        justFirst <- unlist(lapply(strsplit(as.character(myT$Sample.Date.Time), " "), "[[", 1))
        justFirst<-gsub("/16", "/2016", justFirst)
        justFirst <- as.numeric(as.POSIXlt(strptime(justFirst, "%m/%d/%Y")))
        myT$Timestamp <- justFirst

	for( j in (metadataStart + 16):(metadataStart + 25))
	{
		par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
		## bug <- myT[,j]
			justStreams <- myT[(myT$Location == "Mallard Creek" |myT$Location ==  "Sugar Creek") &
						 (myT$Sample == "UP A" |myT$Sample == "UP B" |
						 		 myT$Sample == "DS A" | myT$Sample == "DS B"), ]
                locations <- factor(justStreams$Location)

			streamBugs <- as.numeric(justStreams[,j])
			#myLm <- lm( streamBugs ~ locations)
			names[index] <- names(justStreams)[j]
			pValuesLocation[index] <- wilcox.test( streamBugs[myT$Location == "Mallard Creek"] ,
							 streamBugs[myT$Location == "Sugar Creek"])$p.value
			## mainText = paste(names[index], format(pValuesLocation[index], digits=3))
			boxplot( streamBugs ~ locations, main = paste("Stream (non-parametric) \nuncorrected p-value",  format(pValuesLocation[index], digits=3)))

			updownBinary <- ifelse( justStreams$Sample == "DS A" | justStreams$Sample == "DS B" ,
                                               "down", "up"  )
                updownBinary <- factor(updownBinary, c("up", "down"))

			myFrame <- data.frame(streamBugs,locations, updownBinary, justStreams$Timestamp)

			stripchart(streamBugs ~ locations,
                                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

			fullModel <- lm( streamBugs ~ locations * updownBinary + justStreams$Timestamp )
			anAnova <- anova(fullModel)

			pValueLocationsFromFull[index] <- anAnova$"Pr(>F)"[1]
			pValueUpDownFromFull[index] <- anAnova$"Pr(>F)"[2]
			pValuesTimepointFromFull[index] <- anAnova$"Pr(>F)"[3]
			pValuesUpDownLocationInteraction[index] <- anAnova$"Pr(>F)"[4]

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
                plot(1, type="n", axes=F, xlab="", ylab="")

                        mtext(names[index], outer=TRUE, cex = 1.5)

			index = index + 1
	}

	hist(pValuesLocation,breaks=20, main="Non-parametric Stream")

	hist(pValueLocationsFromFull,breaks=20, main="Stream")
	hist(pValueUpDownFromFull,breaks=20, main="Up/Downstream")
	hist(pValuesTimepointFromFull,breaks=20, main="Timepoint")
	hist(pValuesUpDownLocationInteraction,breaks=20, main="Stream and Up/Down")

	dev.off()
	myFrame <- data.frame(names, pValuesLocation,                              pValueLocationsFromFull,pValueUpDownFromFull,pValuesTimepointFromFull,
	pValuesUpDownLocationInteraction)
	myFrame <- myFrame [order(myFrame$pValuesLocation),]
	myFrame$adjustedPValuesLocation <- p.adjust( myFrame$pValuesLocation, method = "BH" )

	myFrame$adjustedpValueLocationsFromFull <- p.adjust( myFrame$pValueLocationsFromFull, method = "BH" )
	myFrame$adjustedpValueUpDownFromFull <- p.adjust( myFrame$pValueUpDownFromFull, method = "BH" )
	myFrame$adjustedpValuesTimepointFromFull <- p.adjust( myFrame$pValuesTimepointFromFull, method = "BH" )
	myFrame$adjustedpValuesUpDownLocationInteraction <- p.adjust( myFrame$pValuesUpDownLocationInteraction, method = "BH" )

	write.table(myFrame, file=paste("denovoAntibiotic_pValues.txt",sep=""), sep="\t",row.names=FALSE)


