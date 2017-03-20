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

## MDSdata <- read.table("LyteSharon_r01_behavior_pcoa.txt", header=TRUE)

## fullAbundanceMDS
## cecalfecalAbundanceMDS
cecalAbundanceMDS <- read.table("LyteSharon_r01_cecal_pcoa.txt", header=TRUE)

fecalAbundanceMDS <- read.table("LyteSharon_r01_fecal_pcoa.txt", header=TRUE)

fullBehaviorMDS <- read.table("LyteSharon_r01_behavior_pcoa.txt", header=TRUE)

openfieldBehaviorMDS <- read.table("LyteSharon_r01_openfield_pcoa.txt", header=TRUE)
## epmBehaviorMDS
epmBehaviorMDS <- read.table("LyteSharon_r01_epm_pcoa.txt", header=TRUE)
## lightdarkBehaviorMDS
lightdarkBehaviorMDS <- read.table("LyteSharon_r01_lightdark_pcoa.txt", header=TRUE)

MDSlist <- list(fecalAbundanceMDS, fullBehaviorMDS, openfieldBehaviorMDS, epmBehaviorMDS, lightdarkBehaviorMDS)
names(MDSlist) <- c("fecalAbundanceMDS", "fullBehaviorMDS", "openfieldBehaviorMDS", "epmBehaviorMDS", "lightdarkBehaviorMDS")
## May have to go back and do this for cecal contents
myT <- myT[myT$Source == "feces",]
iter <- 0
for(i in MDSlist){
    index <- 1
    iter = iter + 1
    pValuesGroupFromMixed <- vector()
    pValuesSexFromMixed <- vector()
    pValuesSexGroupFromMixedInteraction <- vector()
    pValuesCage <- vector()
    sex <- factor(myT$Sex)
    group <- factor(myT$Exp.or.Ctrl)
    cage <- factor(myT$Cage)
    pdf( paste(names(MDSlist)[iter],".pdf", sep=""))
    for(j in 1:10) {
        par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
        MDSaxis <- i[,j]
        myFrame <- data.frame(MDSaxis, sex, group, cage)

        ## stripchart(MDSaxis ~ sex,
        ##           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

        fullModel <-
            gls( MDSaxis~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

        mixedAnova <- anova(fullModel)
        pValuesGroupFromMixed[index] <- mixedAnova$"p-value"[2]
        pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
        pValuesSexGroupFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

        pValuesCage[index] <- anova(lm( MDSaxis ~ cage ))$"Pr(>F)"[1]
      			boxplot( MDSaxis ~ group,
        main = paste( names(MDSlist)[iter], format(pValuesGroupFromMixed[index],digits=3)))

    stripchart(MDSaxis ~ group ,
               data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names(MDSlist)[iter])
    boxplot( MDSaxis ~ sex,
            main = paste( names(MDSlist)[iter], format(pValuesSexFromMixed[index],digits=3) ) )
    stripchart(MDSaxis ~ sex,
               data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names(MDSlist)[iter])

    boxplot(MDSaxis ~ group*sex, main = paste( names(MDSlist)[iter], format(pValuesSexGroupFromMixedInteraction[index],digits=3) ) )
    stripchart(MDSaxis ~ group*sex,
               data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names(MDSlist)[iter])
    boxplot( MDSaxis ~ cage ,las=2,
            main=format(pValuesCage[index],digits=3))
    stripchart(MDSaxis ~ cage,
               data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names(MDSlist)[iter] )
                               mtext(paste(names(MDSlist)[iter], "axis",j), outer=TRUE, cex = 1.5)

    index = index + 1
}
hist(pValuesGroupFromMixed,breaks=20)
hist(pValuesSexFromMixed,breaks=20)
hist(pValuesSexGroupFromMixedInteraction,breaks=20)
hist(pValuesCage,breaks=20)

dev.off()
## Write out MHC-corrected p-values
myFrame <- data.frame(pValuesGroupFromMixed,
                      pValuesSexFromMixed,pValuesSexGroupFromMixedInteraction, pValuesCage)
myFrame <- myFrame[order(myFrame$pValuesGroupFromMixed),]
myFrame$adjustedpValuesGroupFromMixed<- p.adjust( myFrame$pValuesGroupFromMixed, method = "BH" )
myFrame$adjustedpValuesSexFromMixed<- p.adjust( myFrame$pValuesSexFromMixed, method = "BH" )
myFrame$adjustedpValuesSexGroupFromMixedInteraction<- p.adjust( myFrame$pValuesSexGroupFromMixedInteraction, method = "BH" )
myFrame$pValuesCageAdjust<- p.adjust( myFrame$pValuesCage, method = "BH" )

write.table(myFrame, file=paste(names(MDSlist)[iter], ".txt",sep=""), sep="\t",row.names=FALSE)

}
