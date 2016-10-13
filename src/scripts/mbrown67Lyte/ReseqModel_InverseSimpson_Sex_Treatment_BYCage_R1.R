rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
## Below is what I need for the InverseSimpson diversity.
library("vegan")

## for experiment, treatment, and batch, and no confounder

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")
## All of this data should be present in any analysis, so there is no reason to make it optional or subject to a switch/case statement.

##  inforInverseSimpson<-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE,sep="\t")
## numColsS <- ncol(inforInverseSimpson)
## myColClassesS <- c(rep("character",8), rep("numeric", numColsS-8))
## inforInverseSimpson <-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE,sep="\t",row.names=1,colClasses=myColClassesS)
indexS <- 1
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="")##, na.strings="BLAH")
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
    pValuesInverseSimpson<- vector()
    pValuesKeptCounts<- vector()
    pValuesExperiment <- vector()
    acConf <- vector()
    sexConf <- vector()
    ## index <- 1

    pdf( paste(taxa, "_Reseq_InverseSimpson_Sex_CaseVControl_BYcage_boxplots.pdf", sep=""))

    ## bug <- myT[,i]
    ac <- myT$Condition
    sex <- myT$Sex
    cage <- myT$Cage
    batch <- myT$Sample_Project
    date <- myT$Date
    mo <- myT$MouseOrigin
    treatment <- myT$Treatment
    grouping <- myT$SpreadsheetGrouping
    stressLength <- myT$StressLength
    myT$InverseSimpson <- apply(myT[,2:(ncol(myT) - numMetadataCols)], 1, diversity, index="invsimpson")
    InverseSimpson <- myT$InverseSimpson
    ### myFrame <- data.frame( age, group, sampType, cage, animal, InverseSimpson)
    myFrame <- data.frame(InverseSimpson, ac, sex, cage, treatment, batch, date, mo, grouping, stressLength)

    fullModel <- gls( InverseSimpson ~ sex + treatment, method="REML",correlation=corCompSymm(form=~1|cage),	data = myFrame )
    reducedModel <- gls( InverseSimpson~ sex + treatment, method="REML", data = myFrame )
    fullModelLME <- lme( InverseSimpson~ sex + treatment, method="REML", random = ~1|cage, data = myFrame)

    ## InverseSimpsonPcage[indexS] <- anova(fullModel, reducedModel)$"p-value"[2]
    ## Indicates that the cage effect is not significant in explaining the trend in InverseSimpson diveristy
    ## InverseSimpsonlm <- lm(InverseSimpson ~ group + sampType, x=TRUE)
    ## InverseSimpsonPmodel[indexS] <- list(summary(InverseSimpsonlm)$coefficients[,4])
    ## plot(group, InverseSimpson, main=taxa, ylab="InverseSimpson Diversity")
    ## plot(sampType, InverseSimpson, main=taxa, ylab="InverseSimpson Diversity")

    ## inversesimpson <- myT$inversesimpson
    ## multiWay <- paste(ac, myT$StressLength)
    ## multiWay[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
    ## myFrame <- data.frame(InverseSimpson, ac, sex, cage, treatment, batch, date, mo, grouping, stressLength)

    ## fullModel <- gls( InverseSimpson~  sex + treatment, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
    ## reducedModel <- gls( InverseSimpson~  sex + treatment, method="REML", data = myFrame )
    ## fullModelLME <- lme(InverseSimpson~  sex + treatment, method="REML", random = ~1|factor(cage), data = myFrame)
    ## It seems like reducing the anova calls could speed things up a bit
    ## Remember to revert back to an earlier commit here
    pValuesStressLength[indexS] <- anova(fullModelLME)$"p-value"[5]
    pValuesAcuteChronic[indexS] <- anova(fullModelLME)$"p-value"[5]
    pValuesSex[indexS] <- anova(fullModelLME)$"p-value"[2]
    pValuesExperiment[indexS] <- anova(fullModelLME)$"p-value"[5]
    pValuesBatch[indexS] <- anova(fullModelLME)$"p-value"[5]
    pValuesTreatment[indexS] <- anova(fullModelLME)$"p-value"[3]
    pValuesInteraction[indexS] <- anova(fullModelLME)$"p-value"[5]
    pValuesInverseSimpson[indexS] <- anova(fullModelLME)$"p-value"[5]
    ## Why do you mix model functions here?
    pValuesCage[indexS] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
    iccCage[indexS]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

    ## names[indexS] = names(myT)[i]

    graphMain =  paste( "InverseSimpson Diversity at ", taxa, " levelxs\n",
        ## " pStressLength= ", format(pValuesStressLength[indexS], digits=3), "\n",
        " pSex=", format( pValuesSex[indexS], digits=3),
        ## " pCondition= ", format( pValuesAcuteChronic[indexS],digits=3),
        ## " pBatch= ", format(pValuesBatch[indexS], digits=3), "\n",
        " pTreatment= ", format(pValuesTreatment[indexS], digits=3),"\n",
        ## " pInteraction= ", format(pValuesInteraction[indexS], digits=3),
        ## " pExperiment= ", format(pValuesExperiment[indexS], digits=3),
        ## " pInverseSimpson= ", format(pValuesInverseSimpson[indexS], digits=3),
        ## " pKeptCounts= ", format(pValuesKeptCounts[indexS], digits=3),
        " pCage= " , format( pValuesCage[indexS], digits=3),
        " icc= " , format( iccCage[indexS], digits=3 ), sep="")
    par(mfrow=c(3,1),
        oma = c(1,1,0,0) + 0.1,
        mar = c(1,4,2.5,0) + 0.1)

    plot( InverseSimpson ~ factor(treatment), ylab = "InverseSimpson Diversity", main = graphMain )
    points(factor(sex), InverseSimpson)

    ## plot ( bug ~ factor(ac) )
    ## points(factor(ac), bug)

    plot( InverseSimpson ~ factor(c( paste( myT$Sex, myT$Treatment,sep=""))))
    points(factor(c( paste( myT$Sex, myT$Treatment,sep=""))), InverseSimpson)

    ## plot ( bug ~ factor(treatment) )
    ## points(factor(treatment), bug)

    plot( InverseSimpson ~ factor(cage), ylab=c("InverseSimpson Diversity"))
    points(factor(cage), InverseSimpson)
    indexS=indexS+1

    dFrame <- data.frame( taxa, pValuesSex, pValuesTreatment, pValuesCage, iccCage)#, #pValuesTreatment)#, pValuesBatch)#, pValuesExperiment)
                                        #dropped pValuesSex

    ## dFrame <- dFrame[order(pValuesAcuteChronic),]
    ## dFrame$adjustedStressLength <- p.adjust( dFrame$pValuesStressLength, method = "BH" )
    ## dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
    dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
    ## dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
    dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )
    ## dFrame$adjustedInteraction <- p.adjust( dFrame$pValuesInteraction, method = "BH")
    dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
    ## dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )

    write.table(dFrame,
                file=paste("pValuesForResequencing_InverseSimpson_Sex_Treatment_BYcage_", taxa, ".txt",sep=""),
                sep="\t",row.names=FALSE)

    dev.off()
}
