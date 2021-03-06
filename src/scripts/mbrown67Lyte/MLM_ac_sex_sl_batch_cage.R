rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

# for experiment, treatment, and batch, and no confounder

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")
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
        myColClasses <- c("character", rep("numeric", numCols - numMetadataCols - 1), rep("character", 11), "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character")
       myT <-read.csv(inFileName, header=TRUE, sep="", colClasses=myColClasses, na.strings="BLAH")

#Drop unneeded rows
        removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        myT<-myT[!(myT$Sample_ID %in% removeControls),]
        removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        myT<-myT[!(myT$Sample_ID %in% removetrs),]
        ## removeLow<-c("04-55_S32_L001_R1_001")
        ## myT<-myT[!(myT$MatchFile %in% removeLow),]

                                        #Have to manually drop these for some reason
        manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% manualDrop),]
        labDrop <- c("Harlan Labs")
        myT<-myT[!(myT$MouseOrigin %in% labDrop),]
#print("Got through Processing")
	# our initial model not worrying about confounders except cage
	names <- vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
#       pValuesShannon<- vector()
        pValuesKeptCounts<- vector()
	# WARNING:  EXPERIMENT IS CONFOUNDED WITH AC  + SEX - INTERPRET WITH CAUTION!!!!!
	pValuesExperiment <- vector()
        pValuesStressLength <- vector()
        acConf <- vector()
        sexConf <- vector()
	index <- 1
	pdf( paste(taxa, "_ac_sl_sex_batch_cage_boxplots_DropLab.pdf", sep=""))

	for( i in 2:(ncol(myT) - numMetadataCols))
 		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- myT[,i]
			ac <- myT$Condition
			sex <- myT$Sex
			cage <- myT$Cage
                        batch <- myT$Sample_Project
                        date <- myT$Date
                        mo <- myT$MouseOrigin
                        sl <- myT$StressLength
                        treatment <- myT$Treatment
 #                      shannon <- myT$shannon
                        myFrame <- data.frame(bug, ac, sex, sl, cage, treatment, batch, date, mo)

			fullModel <- gls( bug~   ac + sl + sex + batch, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
                        reducedModel <- gls( bug~  ac + sl + sex + batch, method="REML", data = myFrame )
#                        confCheckModel <- gls (bug~ sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
                        fullModelLME <- lme(bug~  ac + sl + sex + batch, method="REML", random = ~1|factor(cage), data = myFrame)

#                        acConf[index]<-paste((confCheckModel$coefficients["acC"] - fullModel$coefficients["acC"]) / confCheckModel$coefficients["acC"], sep=" ")

#                        sexConf[index]<-paste((confCheckModel$coefficients["sexM"] - fullModel$coefficients["sexM"]) / confCheckModel$coefficients["sexM"], sep=" ")


                                        # It seems like reducing the anova calls could speed things up a bit
                        # Remember to revert back to an earlier commit here
			pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[2]
                        pValuesStressLength[index] <- anova(fullModelLME)$"p-value"[3]
			pValuesSex[index] <- anova(fullModelLME)$"p-value"[4]
			pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[6]
                        pValuesBatch[index] <- anova(fullModelLME)$"p-value"[5]
                        pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[6]
                       # pValuesShannon[index] <- anova(fullModelLME)$"p-value"[5]
                                        #Why do you mix model functions here?
			pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], "\n",
                            " pSex=", format( pValuesSex[index], digits=3),
                            " pCondition= ", format( pValuesAcuteChronic[index],digits=3),
                            " pBatch= ", format(pValuesBatch[index], digits=3), "\n",
                            " pLength= ", format(pValuesStressLength[index], digits=3),
#                            " pTreatment= ", format(pValuesTreatment[index], digits=3),
#                            " pExperiment= ", format(pValuesExperiment[index], digits=3),
                            #" pShannon= ", format(pValuesShannon[index], digits=3),
                            #" pKeptCounts= ", format(pValuesKeptCounts[index], digits=3),
                            " pCage= " , format( pValuesCage[index], digits=3),
                            " icc= " , format( iccCage[index], digits=3 ), sep="")
			par(mfrow=c(4,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)

			plot( bug ~ factor(ac), ylab = names[index],main = graphMain )
                        points(factor(ac), bug)

			plot ( bug ~ factor(batch) )
                       points(factor(batch), bug)

                        plot( bug ~ factor(c( paste( myT$Sex, myT$Condition,sep=""))))
                        points(factor(c( paste( myT$Sex, myT$Condition,sep=""))), bug)

                        plot ( bug ~ factor(sl) )
                        points(factor(sl), bug)

                        plot( bug ~ factor(cage), ylab=names[index])
                        points(factor(cage), bug)
			index=index+1
		}

	dFrame <- data.frame( names, pValuesAcuteChronic, pValuesStressLength, pValuesSex, pValuesBatch, pValuesCage, iccCage)#, #pValuesTreatment)#, pValuesBatch)#, pValuesExperiment)
#dropped pValuesSex
#	dFrame <- dFrame[order(pValuesAcuteChronic),]
	dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
        dFrame$adjustedStressLength <- p.adjust( dFrame$pValuesStressLength, method = "BH" )
	dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
        dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
#       dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )

	dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
#	dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )

        write.table(dFrame, file=paste("pValuesForTaxa_bug_ac_length_sex_batch_cage_DropLab_", taxa, ".txt",sep=""),
                    sep="\t",row.names=FALSE)

        dev.off()
}
