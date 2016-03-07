rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

# for experiment, treatment, and batch, and no confounder

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")
indexS <- 1
ShannonP <- list()
ShannonSummary <- list()

# All of this data should be present in any analysis, so there is no reason to make it optional or subject to a switch/case statement.

#inforShannon<-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE, sep="\t")
#numColsS <- ncol(inforShannon)
#myColClassesS <- c(rep("character",8), rep("numeric", numColsS-8))
#inforShannon <-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE,sep="\t",row.names=1,colClasses=myColClassesS)

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
	myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
	numCols <- ncol(myT)
        numMetadataCols <- 20
        # Reprocessing for correct interpretation
       # myColClasses <- c("character", rep("numeric", numCols - 16), rep("character", 17))
       # myT <-read.csv(inFileName, header=TRUE, sep="", colClasses=myColClasses, na.strings="BLAH")

#Drop unneeded rows
        removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        myT<-myT[!(myT$Sample_ID %in% removeControls),]
        removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        myT<-myT[!(myT$Sample_ID %in% removetrs),]
##        removeLow<-c("04-55_S32_L001_R1_001")
##        myT<-myT[!(myT$MatchFile %in% removeLow),]
        removeDups<-c("04-04_S63_L001_R1_001")
                                        # Somewhat arbitrarily here
        myT<-myT[!(myT$MatchFile %in% removeDups),]

                                        # Have to manually drop these for some reason
        manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% manualDrop),]
        labDrop <- c("Harlan Labs")
        myT<-myT[!(myT$MouseOrigin %in% labDrop),]
	# our initial model not worrying about confounders except cage
        myT$Condition[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
	names <- vector()
        pValuesMultiway<-vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
        pValuesKeptCounts<- vector()
	pValuesExperiment <- vector()
        acConf <- vector()
        multiWay <- vector()
        sexConf <- vector()
	index <- 1
	pdf( paste(taxa, "_testplots.pdf", sep=""))
##         myT$Shannon <- apply(myT[,2:(ncol(myT)-14)],1,diversity)
##         Shannon <- myT$Shannon
##         Day <- myT$Day
##         ImputedBMI<-myT$Imputed.BMI
##         EnergyIntake<-myT$Energy.Intake..kcal.day.

##         Shannonlm <- lm(Shannon ~ Day*patient, x = TRUE)
##         #ShannonlmBMI <- lm(Shannon ~ ImputedBMI*patient, x = TRUE)
##         #ShannonlmEnergy <- lm(Shannon ~ EnergyIntake*patient, x = TRUE)
##         ShannonP[indexS] <- list(summary(Shannonlm)$coefficients[,4][-1])
## #        ShannonSummary[indexS]<-summary(Shannonlm)

## #        names[indexS] = names(myT)[i]

##         indexS <- indexS + 1

##         # Building the data.frames to eventually print out the p-values
##         DayPatientPV.df <- data.frame(ShannonP)
##         DayPatientPV.df <- t(DayPatientPV.df)
##         dFrameDayPatient <- DayPatientPV.df


	for( i in 2:(ncol(myT) - numMetadataCols))
 		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- myT[,i]
			ac <- myT$Condition
			sex <- myT$Sex
			cage <- myT$Cage
                        batch <- myT$Sample_Plate
                        sl <- myT$StressLength
                        date <- myT$Date
                        mo <- myT$MouseOrigin
                        treatment <- myT$Treatment
                        multiWay <- paste(ac, myT$StressLength)
                        multiWay[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"

			names[index] = names(myT)[i]

 			graphMain =  paste( names(myT)[i])
 			par(mfrow=c(3,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)
                        plot(as.numeric(factor(multiWay)), bug, col=as.numeric(batch), pch=16)
                        plot(jitter(as.numeric(sl),.3), bug, col=as.numeric(batch), pch=16 + as.numeric(treatment))
                        plot(as.numeric(batch), bug, col=as.numeric(sex), pch=16)

                        index=index+1
		}
        dev.off()
}
