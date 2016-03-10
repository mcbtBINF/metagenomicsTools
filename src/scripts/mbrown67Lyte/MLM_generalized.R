rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")

for(taxa in taxaLevels){

    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    numMetadataCols <- 20

    ## Reprocessing for correct column data type
    myColClasses <- c("character", rep("numeric", numCols - numMetadataCols), rep("character", 11), "numeric", "numeric", "numeric", "character", "character", "numeric", "character", "character", "character")
    myT <-read.csv(inFileName, header=TRUE, sep="", colClasses=myColClasses, na.strings="BLAH")

    ## Removing unwanted samples {controls, replicates, duplicates, poor quality samples}
    removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
    myT<-myT[!(myT$Sample_ID %in% removeControls),]

    removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
    myT<-myT[!(myT$Sample_ID %in% removetrs),]

    removeDups<-c("04-04_S63_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% removeDups),]

    ## Controls not dropped earlier
    manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% manualDrop),]

    labDrop <- c("Harlan Labs")
    myT<-myT[!(myT$MouseOrigin %in% labDrop),]

    lowSeqDrop <- c("04-55_S32_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% lowSeqDrop),]

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
    pValuesStressLength <- vector()
    pValuesInteraction <- vector()
    acConf <- vector()
    multiWay <- vector()
    allpvals <- list()
    allNames <- list()

    index <- 1

    for( i in 2:(ncol(myT)-numMetadataCols))
        if( sum(myT[,i] != 0 ) > nrow(myT) / 4 ){
            bug <- myT[,i]
            ac <- myT$Condition
            sex <- myT$Sex
            cage <- myT$Cage
            batch <- myT$Sample_Plate
            date <- myT$Date
            mo <- myT$MouseOrigin
            group <- myT$SpreadsheetGrouping
            treatment <- myT$Treatment
            sl <- myT$StressLength
            multiWay <- paste(ac, myT$StressLength)
            multiWay[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
            names[index] = names(myT)[i]

            myFrame <- data.frame(bug, ac, sl, multiWay, sex, cage, treatment, batch, date, mo, date)

            fullModel <- gls( bug~   ac + batch, method="REML",correlation=corCompSymm(form=~1|factor(date)),	data = myFrame )
            reducedModel <- gls( bug~  ac + batch, method="REML", data = myFrame )
            fullModelLME <- lme(bug~  ac + batch, method="REML", random = ~1|factor(date), data = myFrame)
            ## Potential save time by reducing anova calls
            ## Introduce the goodness of fit tests here
            ## Can streamline vectors to be automatically responsive to changes to the model
            allNames[[index]] <- rownames(anova(fullModelLME))[-1]
            allpvals[[index]] <- anova(fullModelLME)$"p-value"[-1]

            pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[2]
            pValuesStressLength[index] <- anova(fullModelLME)$"p-value"[6]
            pValuesInteraction[index] <- anova(fullModelLME)$"p-value"[6]
            pValuesSex[index] <- anova(fullModelLME)$"p-value"[6]
            pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[6]
            pValuesBatch[index] <- anova(fullModelLME)$"p-value"[3]
            pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[6]

            ## Random Effects and Interclass Correlation Coefficient
            pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
            iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

    ## Nice plotting
            ## Correct so that it uses corrected p-values
            ## paste("p-", allNames[[index]], " = ", format(allpvals[[index]], digits=3), " pCage= ", format(pValuesCage[index], digits=3), " icc= " , format( iccCage[index], digits=3 ), sep=""))
            graphMain =  paste( names(myT)[i], "\n",
                ## " pTreatment= ", format(pValuesTreatment[index], digits=3),
                ## " pTime= ", format(pValuesStressLength[index], digits=3),"\n",
                ## " pTreatmentXTime= ", format(pValuesInteraction[index], digits=3),"\n",
                ## " pSex=", format( pValuesSex[index], digits=3),
                " pCondition= ", format( pValuesAcuteChronic[index],digits=3),
                " pBatch= ", format(pValuesBatch[index], digits=3), "\n",
                ## " pTreatment= ", format(pValuesTreatment[index], digits=3),
                ## " pExperiment= ", format(pValuesExperiment[index], digits=3),
                ## " pShannon= ", format(pValuesShannon[index], digits=3),
                ## " pKeptCounts= ", format(pValuesKeptCounts[index], digits=3),
                " pDate= " , format( pValuesCage[index], digits=3),
                " icc= " , format( iccCage[index], digits=3 ), sep="")
            numPlots <- length(allNames[[index]])
            par(mfrow=c(numPlots, 1),
                oma = c(1,1,0,0) + 0.1,
                mar = c(1,4,2.5,0) + 0.1)

            for(i in 1:numPlots){
                if(i == 1){
                    plot( bug ~ factor(get(allNames[[index]][i])), ylab = names[index],main = graphMain )
                    points(factor(get(allNames[[index]][i])), bug)
                }
                else{
                    plot( bug ~ factor(get(allNames[[index]][i])), ylab = names[index])
                    points(factor(get(allNames[[index]][i])), bug)
                }
            }
                ##plot( bug ~ factor(ac), ylab = names[index],main = graphMain )
                ##points(factor(ac), bug)

                ## plot ( bug ~ factor(ac) )
                ## points(factor(ac), bug)

                ## plot( bug ~ factor(c( paste( sl, treatment,sep=""))))
                ## points(factor(c( paste( sl, treatment,sep=""))), bug)

                ##plot ( bug ~ factor(batch) )
                ##points(factor(batch), bug)

                ##plot( bug ~ factor(date), ylab=names[index])
                ##points(factor(date), bug)

            index=index+1
	}

## This should also be dynamic
    dFrame <- data.frame( names, pValuesAcuteChronic, pValuesBatch, pValuesCage, iccCage)
    ## pValuesTreatment)#, pValuesBatch)#, pValuesExperiment)

    ## dFrame <- dFrame[order(pValuesAcuteChronic),]
    dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
    ## dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
    ## dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )
    ## dFrame$adjustedStressLength <- p.adjust( dFrame$pValuesStressLength, method = "BH")
    ## dFrame$adjustedInteraction <- p.adjust( dFrame$pValuesInteraction, method = "BH" )
    dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
    dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
    ## Try and dynamically generate the name of the plot here...
    write.table(dFrame, file=paste("pValuesForTaxa_bug_condition_batch_byDate_", taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    ## Get the sig table working
    ## It would also be nice to have a sig-picker for the plots and a combination of all sig results into one plot.
    keepVector <- grep("adj", names(dFrame))
    sigdFrame <-dFrame[which(dFrame[,keepVector] < 0.05,arr.ind=TRUE),]
    write.table(sigdFrame, file=paste("pValuesForTaxa_bug_condition_batch_byDate_sig_", taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    dev.off()
}
