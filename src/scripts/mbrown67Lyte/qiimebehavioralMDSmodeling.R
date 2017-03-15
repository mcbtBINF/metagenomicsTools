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

MDSdata <- read.table("LyteSharon_r01_behavior_pcoa.txt", header=TRUE)


myT <- myT[myT$Source == "feces",]

	pValuesGroupFromMixed <- vector()
	pValuesSexFromMixed <- vector()
	pValuesSexGroupFromMixedInteraction <- vector()
	pValuesCage <- vector()

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


