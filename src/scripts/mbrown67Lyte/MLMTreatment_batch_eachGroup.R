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
for(groupNum in 2:8){
    for(taxa in taxaLevels ){
	inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
	myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
	numCols <- ncol(myT)
        numMetadataCols <- 20
        ## Reprocessing for correct interpretation
        ## myColClasses <- c("character", rep("numeric", numCols - 16), rep("character", 17))
        ## myT <-read.csv(inFileName, header=TRUE, sep="", colClasses=myColClasses, na.strings="BLAH")
        ## Drop unneeded rows

        removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        myT<-myT[!(myT$Sample_ID %in% removeControls),]
        removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        myT<-myT[!(myT$Sample_ID %in% removetrs),]
        removeDups<-c("04-04_S63_L001_R1_001")
        ## Somewhat arbitrarily here
        myT<-myT[!(myT$MatchFile %in% removeDups),]
        ## removeLow<-c("04-55_S32_L001_R1_001")
        ## myT<-myT[!(myT$MatchFile %in% removeLow),]

        ## Have to manually drop these for some reason
        manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% manualDrop),]

        myT<-myT[myT$SpreadsheetGrouping == groupNum,]
        lowSeqDrop <- c("04-55_S32_L001_R1_001")
        myT<-myT[!(myT$MatchFile %in% lowSeqDrop),]


        ## print("Got through Processing")

	# our initial model not worrying about confounders except cage
	names <- vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
        ## pValuesShannon<- vector()
        pValuesKeptCounts<- vector()
	## WARNING:  EXPERIMENT IS CONFOUNDED WITH AC  + SEX - INTERPRET WITH CAUTION!!!!!
	pValuesExperiment <- vector()
        acConf <- vector()
        sexConf <- vector()
        time <- vector()
        pValuesTime <- vector()

	index <- 1
	pdf( paste(taxa, "_", groupNum,"_Treatment_batch_eachGroup_boxplots.pdf", sep=""))
        print(groupNum)
	for( i in 2:(ncol(myT)-numMetadataCols))
 		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- myT[,i]
			ac <- myT$Condition
			sex <- myT$Sex
			cage <- myT$Cage
                        batch <- myT$Sample_Project
                        date <- myT$Date
                        time <- myT$StressLength
                        mo <- myT$MouseOrigin
                        treatment <- myT$Treatment

                        myFrame <- data.frame(bug, ac, sex, cage, treatment, batch, date, mo)

                        fullModel <- lm(bug ~ treatment + batch, x=TRUE)

			##fullModel <- gls( bug~   treatment + batch, method="REML",correlation=corCompSymm(form=~1|factor(cage)), data = myFrame )
                        ##reducedModel <- gls( bug~  treatment + batch, method="REML", data = myFrame )
                        ##fullModelLME <- lme(bug~  treatment + batch, method="REML", random = ~1|factor(cage), data = myFrame)
                        ## It seems like reducing the anova calls could speed things up a bit
                        pValuesTreatment[index] <- anova(fullModel)$"Pr(>F)"[1]
                        pValuesBatch[index] <- anova(fullModel)$"Pr(>F)"[2]
      			##pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[2]
                        ##pValuesTime[index] <- anova(fullModelLME)$"p-value"[5]
			## pValuesSex[index] <- anova(fullModelLME)$"p-value"[5]
			## pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[5]
                        ## pValuesStressLength[index] <- anova(fullModel)$"Pr(>F)"[2]
                        ##pValuesBatch[index] <- anova(fullModelLME)$"p-value"[3]

                        ## Why do you mix model functions here?

			## pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			## iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], "\n",
                            ## " pSex=", format( pValuesSex[index], digits=3),
                            ## " pAcuteChronic= ", format( pValuesAcuteChronic[index],digits=3),
                            " pBatch= ", format(pValuesBatch[index], digits=3),
                            " pTreatment= ", format(pValuesTreatment[index], digits=3)#, "\n",
                            ## " pExperiment= ", format(pValuesExperiment[index], digits=3),
                            ##" pCage= " , format( pValuesCage[index], digits=3),
                            ##" icc= " , format( iccCage[index], digits=3 ), sep=""
                                           )
			par(mfrow=c(2,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)

			plot( bug ~ factor(batch), ylab = names[index],main = graphMain )
                        points(factor(batch), bug)

			## plot ( bug ~ factor(ac) )
                        ## points(factor(ac), bug)

                        plot( bug ~ factor(c( paste( myT$Sex, treatment,sep=""))))
                        points(factor(c( paste( myT$Sex, treatment,sep=""))), bug)

                        ## plot ( bug ~ factor(treatment) )
                        ## points(factor(treatment), bug)

                        ## plot( bug ~ factor(cage), ylab=names[index])
                        ## points(factor(cage), bug)
			index = index+1
		}

	dFrame <- data.frame( names, pValuesTreatment, pValuesBatch)#, pValuesCage, iccCage)#, #pValuesTreatment)#, pValuesBatch)#, pValuesExperiment)

	## dFrame <- dFrame[order(pValuesAcuteChronic),]
	## dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
	dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )
        dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
        ## dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )

##	dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
	## dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )

        write.table(dFrame, file=paste("pValuesForTaxa_bug_treatment_batch_eachgroup_", taxa, "_", groupNum, ".txt",sep=""),
                    sep="\t",row.names=FALSE)

##        keepVector <- grep("adj", names(dFrame))
##        sigdFrame <-dFrame[which(dFrame[,keepVector] < 0.05,arr.ind=TRUE),]
##        write.table(sigdFrame, file=paste("pValuesForTaxa_bug_treatment_batch_nocage_eachgroup_sig_", taxa, "_", groupNum,".txt",sep=""), sep="\t",row.names=FALSE)

        dev.off()
    }
}
