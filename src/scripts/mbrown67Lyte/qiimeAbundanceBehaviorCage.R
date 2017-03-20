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

myT <- myT[myT$Source == "feces",]
iter <- 1
    pValuesGroupFromMixed <- vector()
    pValuesSexFromMixed <- vector()
    pValuesSexGroupFromMixedInteraction <- vector()
    pValuesCage <- vector()
names <- list()
    pdf( paste("BehavioralAssociation",".pdf", sep=""))
index <- 1
	for ( j in c(endMetadataIndex + 1, ncol(myT)))	{
		## par(mfrow=c(2,2))
		bug <- myT[,j]

		if( sum(bug != 0 ) > nrow(subT) / 4 )
		{

for(i in 20:58){
    ## iter = iter + 1
    sex <- factor(myT$Sex)
    group <- factor(myT$Exp.or.Ctrl)
    cage <- factor(myT$Cage)
## par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
names[index]<-colnames(myT)[i]
behave <- myT[,i]
myFrame <- data.frame(behave, sex, group, cage)
        fullModel <-
            gls( bug ~ behave, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

        mixedAnova <- anova(fullModel)
        pValuesGroupFromMixed[index] <- mixedAnova$"p-value"[2]
        ## pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
        ## pValuesSexGroupFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

        pValuesCage[index] <- anova(lm( behave ~ cage ))$"Pr(>F)"[1]
    ##   			boxplot( behave ~ group,
    ##     main = paste( names[index], format(pValuesGroupFromMixed[index],digits=3)))

    ## stripchart(behave ~ group ,
    ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
    ## boxplot( behave ~ sex,
    ##         main = paste( names[index], format(pValuesSexFromMixed[index],digits=3) ) )
    ## stripchart(behave ~ sex,
    ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

    ## boxplot(behave ~ group*sex, main = paste( names[index], format(pValuesSexGroupFromMixedInteraction[index],digits=3) ) )
    ## stripchart(behave ~ group*sex,
    ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
    ## boxplot( behave ~ cage ,las=2,
    ##         main=format(pValuesCage[index],digits=3))
    ## stripchart(behave ~ cage,
    ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )
    ##                            mtext(names[index], outer=TRUE, cex = 1.5)
index = index + 1
}

hist(pValuesGroupFromMixed,breaks=20)
hist(pValuesSexFromMixed,breaks=20)
hist(pValuesSexGroupFromMixedInteraction,breaks=20)
hist(pValuesCage,breaks=20)

dev.off()
## Write out MHC-corrected p-values
myFrame <- data.frame(colnames(myT)[20:58], pValuesGroupFromMixed,
                      pValuesSexFromMixed,pValuesSexGroupFromMixedInteraction, pValuesCage)
myFrame <- myFrame[order(myFrame$pValuesGroupFromMixed),]
myFrame$adjustedpValuesGroupFromMixed<- p.adjust( myFrame$pValuesGroupFromMixed, method = "BH" )
myFrame$adjustedpValuesSexFromMixed<- p.adjust( myFrame$pValuesSexFromMixed, method = "BH" )
myFrame$adjustedpValuesSexGroupFromMixedInteraction<- p.adjust( myFrame$pValuesSexGroupFromMixedInteraction, method = "BH" )
myFrame$pValuesCageAdjust<- p.adjust( myFrame$pValuesCage, method = "BH" )

write.table(myFrame, file=paste("AbundanceBehaviorCagepValues", ".txt",sep=""), sep="\t",row.names=FALSE)
