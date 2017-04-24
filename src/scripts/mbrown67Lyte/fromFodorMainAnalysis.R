rm(list=ls())

## setwd("C:\\LyteBehaviorMarch2017\\rg_results\\")
library("Kendall")
library("pscl")
library("lmtest")
library("nlme")

## myT <- read.table("LyteSharon_r01_crDataOnlyTaxaAsColumnsLogNormPlusMetadata.txt", sep="\t", header=TRUE)

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

inFileName <- "qiimeOTULogNormwithMetadata.txt"
myT <- read.table(inFileName,header=TRUE,sep="\t")

otuMapping <- read.table("LyteSharon_r01_cr_MAPPING.txt", sep="\t", header=TRUE)
otuMapping[,1] <- paste0("X",otuMapping[,1])
rownames(otuMapping)<-otuMapping[,1]

## tissues <- unique(myT$Source)
tissues <- c("feces")

for( t in tissues ) {
    subT <- myT[ myT$Source == t, ]
    index <- 1
    names <- vector()
    pValuesGroup <- vector()
    pValuesSex <- vector()

    pValuesGroupFromMixed <- vector()
    pValuesSexFromMixed <- vector()
    pValuesSexGroupFromMixedInteraction <- vector()
    pValuesCage <- vector()

    ## pdf( paste("qiimeOTUanalysis_", t,".pdf", sep=""))
    pdf( paste("qiimeOTUanalysis_", t,".pdf", sep=""))


    ## Ignores the metadata
    for ( j in c(61:3560)) {
        ## par(mfrow=c(2,2), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
        bug <- subT[,j]

        if( sum(bug != 0 ) > nrow(subT) / 4 ) {
            pValuesGroup[index] <- wilcox.test( subT[subT$Group=="Experimental", j],
                                               subT[subT$Group=="Control", j])$p.value
            names[index] <- names(subT)[j]

            group <- factor(subT$Group)
            sex <- factor(subT$Sex)

            cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                       "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))

            myFrame <- data.frame(bug, sex, group, cage)

            pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
                                             subT[subT$Sex=="female", j])$p.value

            fullModel <-
                gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

            mixedAnova <- anova(fullModel)
            pValuesGroupFromMixed[index] <- mixedAnova$"p-value"[2]
            pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
            pValuesSexGroupFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

            pValuesCage[index] <- anova(lm( bug ~ cage ))$"Pr(>F)"[1]

            index <- index + 1
        }

    }
    index <- 1

    namesMapping <- otuMapping[names,2]

    dFrame <- data.frame(names,namesMapping, pValuesGroup,pValuesSex,pValuesGroupFromMixed,
                         pValuesSexFromMixed,pValuesSexGroupFromMixedInteraction, pValuesCage)
    dFrame <- dFrame [order(dFrame$pValuesGroup),]
    dFrame$adjustedPValuesGroup <- p.adjust( dFrame$pValuesGroup, method = "BH" )
    dFrame$adjustedPValuesSex<- p.adjust( dFrame$pValuesSex, method = "BH" )

    dFrame$adjustedpValuesGroupFromMixed<- p.adjust( dFrame$pValuesGroupFromMixed, method = "BH" )
    dFrame$adjustedpValuesSexFromMixed<- p.adjust( dFrame$pValuesSexFromMixed, method = "BH" )
    dFrame$adjustedpValuesSexGroupFromMixedInteraction<- p.adjust( dFrame$pValuesSexGroupFromMixedInteraction, method = "BH" )
    dFrame$pValuesCageAdjust<- p.adjust( dFrame$pValuesCage, method = "BH" )
    rownames(dFrame) <- dFrame$names

    write.table(dFrame, file=paste("qiimeOTUanalysis", t, ".txt",sep=""), sep="\t",row.names=FALSE)

    ## for ( j in c(61:3560)) {
    for (nameIter in rownames(dFrame)) {
        ## This will need to be tweaked
        ## par(mfrow=c(2,2), mgp=c(3, 1.25, 0)) ##, oma = c(0, 0, 2, 0))
        par(mfrow=c(2,2),
            oma = c(0, 3, 2, 0) + 0.1,
            mar = c(3, 2, 2, 1) + 0.1,
            mgp = c(3, 1.25, 0))

        bug <- subT[,nameIter]

        ## if( sum(bug != 0 ) > nrow(subT) / 4 ) {
        ## names[index] <- names(subT)[j]

        group <- factor(subT$Group)
        sex <- factor(subT$Sex)

        cage <-  factor( paste( subT$Housing,s ubT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                   "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))

        myFrame <- data.frame(bug, sex, group, cage)

        boxplot( bug ~ group , ylab="Log normalized abundance",
                main = paste("stress p-value", format(dFrame$adjustedpValuesGroupFromMixed[index],digits=3) ) )

        stripchart(bug ~ group ,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = "Log normalized abundance")

        boxplot( bug ~ sex,
                main = paste("sex p-value", format(dFrame$adjustedpValuesSexFromMixed[index],digits=3) ), ylab="Log normalized abundance" )

        stripchart(bug ~ sex,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

        boxplot( bug ~ group*sex, ylab="Log normalized abundance",
                main=paste("stress*sex p-value", format(dFrame$adjustedpValuesSexGroupFromMixedInteraction[index],digits=3)),
                xaxt='n')
        axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)

        stripchart(bug ~ group*sex,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )

        boxplot( bug ~ cage, ylab="Log normalized abundance",
                main=paste("cage effect p-value", format(dFrame$pValuesCageAdjust[index],digits=3)), xaxt='n')
        ##            axis(1, at=c(1, 2, 3, 4, 5, 6, 7, 8),
        ##                 labels=c("FC1", "FC2", "FX1", "FX2", "MC1", "MC2", "MX1", "MX2"),
        ##                 cex.axis=0.65)
        axis(1, at=c(1, 3, 5, 7),
             labels=c("FC1", "FX1", "MC1", "MX1"), cex.axis=0.8)
        axis(1, at=c(2, 4, 6, 8),
             labels=c("FC2", "FX2", "MC2", "MX2"), cex.axis=0.8)

        stripchart(bug ~ cage,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])
        ## Fix this for the presentation, but consider changing it back later
        mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)
        mtext("Log normalized abundance", outer=TRUE, side=2, line=1)
        ## plot(1, type="n", axes=F, xlab="", ylab="")

        index = index + 1
    }

    hist(pValuesGroupFromMixed,breaks=20, main="")
    hist(pValuesSexFromMixed,breaks=20, main="")
    hist(pValuesSexGroupFromMixedInteraction,breaks=20, main="")
    hist(pValuesCage,breaks=20, main="")

    dev.off()
    ## Code for plotting volcano plots (4) here
    ## Need the p-values
    ## Need the average abundance before and
}
