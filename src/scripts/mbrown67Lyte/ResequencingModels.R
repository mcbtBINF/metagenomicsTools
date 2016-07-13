rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

# for experiment, treatment, and batch, and no confounder

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")
# All of this data should be present in any analysis, so there is no reason to make it optional or subject to a switch/case statement.

#inforShannon<-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE, sep="\t")
#numColsS <- ncol(inforShannon)
#myColClassesS <- c(rep("character",8), rep("numeric", numColsS-8))
#inforShannon <-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE,sep="\t",row.names=1,colClasses=myColClassesS)

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
	myT <-read.csv(inFileName,header=TRUE,sep="")#, na.strings="BLAH")
	numCols <- ncol(myT)
        numMetadataCols <- 19

        ## Dropping Pos/Neg and samples from mice from the other lab
        myT <- myT[-which(is.na(myT$MouseOrigin) == TRUE),]
        myT <- myT[-which(myT$MouseOrigin == "Harlan Labs"), ]

	## our initial model not worrying about confounders except cage
        ## myT$Condition[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
	names <- vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        pValuesStressLength <- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
        pValuesInteraction <- vector()
        ##       pValuesShannon<- vector()
        pValuesKeptCounts<- vector()
	pValuesExperiment <- vector()
        acConf <- vector()
        sexConf <- vector()
	index <- 1

	pdf( paste(taxa, "_ReseqTreatment_Sex_BYcage_boxplots.pdf", sep=""))

	for( i in 2:(ncol(myT)-numMetadataCols))
 		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- myT[,i]
			ac <- myT$Condition
			sex <- myT$Sex
			cage <- myT$Cage
                        batch <- myT$Sample_Project
                        date <- myT$Date
                        mo <- myT$MouseOrigin
                        treatment <- myT$Treatment
                        grouping <- myT$SpreadsheetGrouping
                        stressLength <- myT$StressLength
                        ## shannon <- myT$shannon
                        ## multiWay <- paste(ac, myT$StressLength)
                        ## multiWay[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
                        myFrame <- data.frame(bug, ac, sex, cage, treatment, batch, date, mo, grouping, stressLength)

			fullModel <- gls( bug~  treatment + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
                        reducedModel <- gls( bug~  treatment + sex, method="REML", data = myFrame )
                        fullModelLME <- lme(bug~  treatment + sex, method="REML", random = ~1|factor(cage), data = myFrame)
                                        # It seems like reducing the anova calls could speed things up a bit
                                        # Remember to revert back to an earlier commit here
                        pValuesStressLength[index] <- anova(fullModelLME)$"p-value"[5]
			pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[5]
			pValuesSex[index] <- anova(fullModelLME)$"p-value"[3]
			pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[5]
                        pValuesBatch[index] <- anova(fullModelLME)$"p-value"[5]
                        pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[2]
                        pValuesInteraction[index] <- anova(fullModelLME)$"p-value"[5]
                       # pValuesShannon[index] <- anova(fullModelLME)$"p-value"[5]
                                        #Why do you mix model functions here?
			pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], "\n",
##                            " pStressLength= ", format(pValuesStressLength[index], digits=3), "\n",
                            " pSex=", format( pValuesSex[index], digits=3),
#                            " pCondition= ", format( pValuesAcuteChronic[index],digits=3),
#                            " pBatch= ", format(pValuesBatch[index], digits=3), "\n",
                            " pTreatment= ", format(pValuesTreatment[index], digits=3),"\n",
#                            " pInteraction= ", format(pValuesInteraction[index], digits=3),
#                            " pExperiment= ", format(pValuesExperiment[index], digits=3),
                            #" pShannon= ", format(pValuesShannon[index], digits=3),
                            #" pKeptCounts= ", format(pValuesKeptCounts[index], digits=3),
                            " pCage= " , format( pValuesCage[index], digits=3),
                            " icc= " , format( iccCage[index], digits=3 ), sep="")
			par(mfrow=c(3,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)

			plot( bug ~ factor(treatment), ylab = names[index],main = graphMain )
                        points(factor(sex), bug)

#			plot ( bug ~ factor(ac) )
#                       points(factor(ac), bug)

                        plot( bug ~ factor(c( paste( myT$Sex, myT$Treatment,sep=""))))
                        points(factor(c( paste( myT$Sex, myT$Treatment,sep=""))), bug)

#                        plot ( bug ~ factor(treatment) )
#                        points(factor(treatment), bug)

                        plot( bug ~ factor(cage), ylab=names[index])
                        points(factor(cage), bug)
			index=index+1
		}

	dFrame <- data.frame( names, pValuesTreatment, pValuesSex, pValuesCage, iccCage)#, #pValuesTreatment)#, pValuesBatch)#, pValuesExperiment)
                                        #dropped pValuesSex

#	dFrame <- dFrame[order(pValuesAcuteChronic),]
##        dFrame$adjustedStressLength <- p.adjust( dFrame$pValuesStressLength, method = "BH" )
                                        #	dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
	dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
#        dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
        dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )
#        dFrame$adjustedInteraction <- p.adjust( dFrame$pValuesInteraction, method = "BH")
	dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
#	dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )

        write.table(dFrame, file=paste("pValuesForResequencing_bug_Treatment_Sex_BYcage_", taxa, ".txt",sep=""),
                    sep="\t",row.names=FALSE)

        dev.off()
}
