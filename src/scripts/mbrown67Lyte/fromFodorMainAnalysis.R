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

tissues <- unique(myT$Source)

for( t in tissues )
{
	subT <- myT[ myT$Source == t, ]
	index <- 1
	names <- vector()
	pValuesGroup <- vector()
	pValuesSex <- vector()

	pValuesGroupFromMixed <- vector()
	pValuesSexFromMixed <- vector()
	pValuesSexGroupFromMixedInteraction <- vector()
	pValuesCage <- vector()

	pdf( paste("qiimeOTUanalysis_", t,".pdf", sep=""))

        ## Ignores the metadata
        for ( j in c(61:3560))
	{
		par(mfrow=c(2,2))
		bug <- subT[,j]

		if( sum(bug != 0 ) > nrow(subT) / 4 )
		{
			pValuesGroup[index] <- wilcox.test( subT[subT$Group=="Experimental", j],
										subT[subT$Group=="Control", j])$p.value
			names[index] <- names(subT)[j]

			group <- factor(subT$Group)
			boxplot( bug ~ group ,
			main = paste( names[index], format(pValuesGroup[index],digits=3) ) )

			myFrame <- data.frame(bug, group)

			stripchart(bug ~ group ,
				data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

			sex <- factor(subT$Sex);

			pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
										subT[subT$Sex=="female", j])$p.value

			boxplot( bug ~ sex,
 			main = paste( names[index], format(pValuesSex[index],digits=3) ) )

			cage <-  factor( paste( subT$Housing, subT$Sex, sep=""))


			myFrame <- data.frame(bug, sex,group,cage)

			stripchart(bug ~ sex,
				data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

			fullModel <-
			gls( bug~  group * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

			mixedAnova <- anova(fullModel)
			pValuesGroupFromMixed[index] <- mixedAnova$"p-value"[2]
			pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
			pValuesSexGroupFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

			pValuesCage[index] <- anova(lm( bug ~ cage ))$"Pr(>F)"[1]
                        boxplot( bug ~ group*sex ,las=2,
				main=format(pValuesSexGroupFromMixedInteraction[index],digits=3))
			stripchart(bug ~ group*sex,
				data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )

			boxplot( bug ~ cage ,las=2,
				main=format(pValuesCage[index],digits=3))
			stripchart(bug ~ cage,
				data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )

			## plot(1, type="n", axes=F, xlab="", ylab="")

			index <- index + 1
		}

	}

	## hist(pValuesGroup,breaks=20)
	## hist(pValuesSex,breaks=20)
	hist(pValuesGroupFromMixed,breaks=20)
	hist(pValuesSexFromMixed,breaks=20)
	hist(pValuesSexGroupFromMixedInteraction,breaks=20)
	hist(pValuesCage,breaks=20)

	dFrame <- data.frame(names,pValuesGroup,pValuesSex,pValuesGroupFromMixed,
	pValuesSexFromMixed,pValuesSexGroupFromMixedInteraction, pValuesCage)
	dFrame <- dFrame [order(dFrame$pValuesGroup),]
	dFrame$adjustedPValuesGroup <- p.adjust( dFrame$pValuesGroup, method = "BH" )
	dFrame$adjustedPValuesSex<- p.adjust( dFrame$pValuesSex, method = "BH" )

	dFrame$adjustedpValuesGroupFromMixed<- p.adjust( dFrame$pValuesGroupFromMixed, method = "BH" )
	dFrame$adjustedpValuesSexFromMixed<- p.adjust( dFrame$pValuesSexFromMixed, method = "BH" )
	dFrame$adjustedpValuesSexGroupFromMixedInteraction<- p.adjust( dFrame$pValuesSexGroupFromMixedInteraction, method = "BH" )
	dFrame$pValuesCageAdjust<- p.adjust( dFrame$pValuesCage, method = "BH" )

	write.table(dFrame, file=paste("qiimeOTUanalysis", t, ".txt",sep=""), sep="\t",row.names=FALSE)
	dev.off()


}
