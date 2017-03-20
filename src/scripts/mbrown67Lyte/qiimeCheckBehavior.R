rm(list=ls())

## Load modeling packages
library("Kendall")
library("pscl")
library("lmtest")
library("nlme")

## load data
## load MDS data

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

inFileName <- "qiimeOTULogNormwithMetadata.txt"
myT <- read.table(inFileName,header=TRUE,sep="\t")

endMetadataIndex <- which(colnames(myT) == "counts")

    index <- 1
    pValuesGroupFromMixed <- vector()
    pValuesSexFromMixed <- vector()
    pValuesSexGroupFromMixedInteraction <- vector()
    pValuesCage <- vector()
    sex <- myT$Sex
    group <- myT$Exp.or.Ctrl
    cage <- myT$Cage

for(behavior in 20:58){
    behave <- myT[,behavior]
    myFrame <- data.frame(behave, sex, group, cage)
            fullModel <-
            gls( behave~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

        mixedAnova <- anova(fullModel)
        pValuesGroupFromMixed[index] <- mixedAnova$"p-value"[2]
        pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
        pValuesSexGroupFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

        pValuesCage[index] <- anova(lm( behave ~ cage ))$"Pr(>F)"[1]
    index = index + 1
    }
        	myFrame <- data.frame(pValuesGroupFromMixed,
	pValuesSexFromMixed,pValuesSexGroupFromMixedInteraction, pValuesCage)
	myFrame <- myFrame[order(myFrame$pValuesGroupFromMixed),]
	myFrame$adjustedpValuesGroupFromMixed<- p.adjust( myFrame$pValuesGroupFromMixed, method = "BH" )
	myFrame$adjustedpValuesSexFromMixed<- p.adjust( myFrame$pValuesSexFromMixed, method = "BH" )
	myFrame$pValuesSexGroupFromMixedInteraction<- p.adjust( myFrame$pValuesSexGroupFromMixedInteraction, method = "BH" )
	myFrame$pValuesCageAdjust<- p.adjust( myFrame$pValuesCage, method = "BH" )

	## write.table(myFrame, file=paste("CheckBehavior", ".txt",sep=""), sep="\t",row.names=FALSE)


