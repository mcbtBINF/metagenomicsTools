rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
	myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
	numCols <- ncol(myT)
        numMetadataCols <- 20

        ## Dropping unneeded samples
        removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        myT<-myT[!(myT$Sample_ID %in% removeControls),]
        removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        myT<-myT[!(myT$Sample_ID %in% removetrs),]
        ## Remove sample duplicated across both batches
        removeDups<-c("04-04_S63_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% removeDups),]
        ## Have to manually drop these for some reason
        manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% manualDrop),]
        ## Dropping the protype experiment as it is from a different mouse source.
        ## Other analyses support this decision.
        labDrop <- c("Harlan Labs")
        myT<-myT[!(myT$MouseOrigin %in% labDrop),]

        ## Renaming so that control mice are redundantly labeled.
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

                        myFrame <- data.frame(bug, ac, multiWay, sex, cage, treatment, batch, date, mo)

			## fullModel <- gls( bug~   ac, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
                        ## vals<-list(summary(fullModel)$coefficients[,4][-1])

			names[index] = names(myT)[i]

 			graphMain =  paste( names(myT)[i])
 			par(mfrow=c(5,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)

      			fullModel1 <- gls( bug~ sl*treatment, method="REML",correlation=corCompSymm(form=~1|factor(cage)), data = myFrame )

                        plot(jitter(as.numeric(sl), .3), bug, col=as.numeric(batch) + 3, pch=16 + as.numeric(treatment), xlab= "14, 19, 9")

                        abline(a = fullModel1$coef[1], b = fullModel1$coef[2])
                        abline(a = fullModel1$coef[1] + fullModel1$coef[3], b = fullModel1$coef[2] + fullModel1$coef[5], col="RED")
                        abline(a = fullModel1$coef[1] + fullModel1$coef[4], b = fullModel1$coef[2] + fullModel1$coe[6], col="GREEN")
                        ##abline(a = fullModel1$coef[1] + fullModel1$coef[4], b = fullModel1$coef[6] + fullModel1$coef[2], col="BLUE")

                        plot(jitter(as.numeric(factor(multiWay)),.3), bug, col=as.numeric(factor(multiWay)), pch=16 + as.numeric(batch), xlab= "Acute_9, Chronic_14, Chronic_19, Control")


                        plot(jitter(as.numeric(sl), .3), bug, col=as.numeric(batch), pch=16 + as.numeric(treatment), xlab= "14, 19, 9")
                        plot(jitter(as.numeric(batch).3), bug, col=as.numeric(sex), pch=16, xlab="Which batch?")
                        plot(jitter(as.numeric(myT$SpreadsheetGrouping),.3), bug, xlab="Spreadsheet Group number", pch = 16)
                        plot(jitter(as.numeric(treatment),.3), bug, xlab = "Control versus Case", pch=16)
                        plot(jitter(as.numeric(date),.3), bug, pch=16, xlab="11/30/14, 2/17/15, 3/18/15, 4/14/15")
                        ## axis(1, at=as.numeric(date), labels = c("11/30/14", "2/17/15", "3/18/15", "4/14/15"))


                        index=index+1

                        ## Try for a plot with abline from a simple linear fit
                    }

        dev.off()
}
