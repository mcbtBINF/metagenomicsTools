rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

# for keptCounts with experiment

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/FodorProcessedForwardReadsNoCutoff/noLowTaxas")
#setwd("C:\\MarkRatDataAug2015")

taxaLevels <- c( "p", "c", "o", "f", "g", "otu");
# All of this data should be present in any analysis, so there is no reason to make it optional or subject to a switch/case statement.
raw.Freads<-read.table("../../qiime18_Freads_cr.txt", header=TRUE, comment.char="@", sep="\t")
colcounts<-colSums(raw.Freads[,2:70])
#boxplot(colcounts)
#points(rnorm(69,mean=1,sd=0.05), colcounts)
# Drop low reads (less than 1000) and negative controls and water and technical replicates
excludedCols<-c("X.OTUID", "X45", "X50", "X59tr", "X60tr", "X80", "X81", "X82", "X83", "Water", "N1", "N2", "N3", "P2", "taxonomy")
keptSamples<-raw.Freads[ , -which(names(raw.Freads) %in% excludedCols)]
keptCounts<-colSums(keptSamples)

inforShannon<-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE, sep="\t")
numColsS <- ncol(inforShannon)
myColClassesS <- c(rep("character",8), rep("numeric", numColsS-8))
inforShannon <-read.table("otuTaxaAsColumnsLogNormWithMetadata.txt.temp", header=TRUE,sep="\t",row.names=1,colClasses=myColClassesS)

for(taxa in taxaLevels )
{
	inFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetadata.txt.tempC", sep ="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c(rep("character",8), rep("numeric", numCols-8))
	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)

	# squash the negative controls
	myT <- myT[ which( myT$sex != "Non"), ]

# get rid of low depth samples
#        myT$shannon<-diversity(inforShannon[,8:(numColsS - 1)], index="shannon", MARGIN=1)
        myT$keptCounts<-keptCounts

	# our initial model not worrying about confounders except cage
	names <- vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
#        pValuesShannon<- vector()
        pValuesKeptCounts<- vector()
	# WARNING:  EXPERIMENT IS CONFOUNDED WITH AC  + SEX - INTERPRET WITH CAUTION!!!!!
	pValuesExperiment <- vector()
	index <- 1
	pdf( paste(taxa, "boxplots.pdf", sep=""))

	for( i in 9:ncol(myT))
            #This is the poor OTU/resolution limiting step
		if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
		{
			bug <- myT[,i]
			ac <- myT$acuteOrChronic
			sex <- myT$sex
			cage <- myT$cage
			experiment <- myT$expriment
                        batch <- myT$batch
                        treatment <- myT$treatment
#                        shannon <- myT$shannon
                        keptCounts <- myT$keptCounts
                        myFrame <- data.frame(keptCounts, ac, sex, cage, experiment, batch, treatment)
#                        print(dim(myFrame))
                                        #myFrame <- data.frame(shannon, ac, sex, cage, experiment, batch, treatment)#,
                        #myFrame <- data.frame(bug, ac, sex, cage, experiment, batch, treatment)#, keptCounts)#, shannon, keptCounts)
#                        myFrame <- data.frame(keptCounts, ac, sex, cage, experiment, batch, treatment)#, keptCounts)#, shannon, keptCounts)
#			myFrame <- data.frame(bug, ac, sex, cage, experiment, batch, treatment)#, keptCounts)#, shannon, keptCounts)
#			myFrame <- data.frame(bug, ac, sex, cage, experiment, batch, treatment)#, keptCounts
#              		fullModel <- gls( shannon~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
#			fullModel <- gls( keptCounts~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
#			fullModel <- gls( bug~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
#			fullModel <- gls( bug~  ac + sex + experiment, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
			fullModel <- gls( keptCounts~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
#			fullModel <- gls( bug~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(batch)),				data = myFrame )
#			reducedModel <- gls( shannon~  ac + sex, method="REML",	data = myFrame )
#			reducedModel <- gls( keptCounts~  ac + sex, method="REML",	data = myFrame )
                                        #			reducedModel <- gls( bug~  ac + sex, method="REML",	data = myFrame )
                        reducedModel <- gls( keptCounts~  ac + sex, method="REML",	data = myFrame )
#			reducedModel <- gls( bug~  ac + sex + experiment, method="REML",	data = myFrame )
#			reducedModel <- gls( bug~  ac + sex, method="REML",	data = myFrame )

#			fullModelLME <- lme(shannon~  ac + sex, method="REML", random = ~1|factor(cage), data = myFrame)
#			fullModelLME <- lme(keptCounts~  ac + sex, method="REML", random = ~1|factor(cage), data = myFrame)
#			fullModelLME <- lme(bug~  ac + sex, method="REML", random = ~1|factor(cage), data = myFrame)
#			fullModelLME <- lme(bug~  ac + sex + experiment, method="REML", random = ~1|factor(cage), data = myFrame)
                        fullModelLME <- lme(keptCounts~  ac + sex, method="REML", random = ~1|factor(cage), data = myFrame)
#			fullModelLME <- lme(bug~  ac + sex, method="REML", random = ~1|factor(batch), data = myFrame)

# It seems like reducing the anova calls could speed things up a bit
			pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[2]
			pValuesSex[index] <- anova(fullModelLME)$"p-value"[3]
			pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[5]
                        pValuesBatch[index] <- anova(fullModelLME)$"p-value"[5]
                        pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[5]
 #                       pValuesShannon[index] <- anova(fullModelLME)$"p-value"[5]
                                        #Why do you mix model functions here?
			pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
			iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

			names[index] = names(myT)[i]

			graphMain =  paste( names(myT)[i], "\n",
                            " pSex=", format( pValuesSex[index], digits=3),
                            " pAcuteChronic= ", format( pValuesAcuteChronic[index],digits=3),
#                            " pBatch= ", format(pValuesBatch[index], digits=3), "\n",
#                            " pTreatment= ", format(pValuesTreatment[index], digits=3),
#                            " pExperiment= ", format(pValuesExperiment[index], digits=3),
                            #" pShannon= ", format(pValuesShannon[index], digits=3),
                            #" pKeptCounts= ", format(pValuesKeptCounts[index], digits=3),
                            " pCage= " , format( pValuesCage[index], digits=3),
                            " icc= " , format( iccCage[index], digits=3 ), sep="")
                        #layout(matrix(1:5, ncol = 1), widths = 1,
                        #       heights = c(2,2,2,2,2.5), respect = FALSE)
			par(mfrow=c(4,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)
                        #figure margins too large...

#                       plot( shannon ~ factor(sex), ylab = names[index],main = graphMain )
#                        plot( keptCounts ~ factor(sex), ylab = names[index],main = graphMain )
                                        #                        par(mar=c(1,3,3,2)+0.1)
#                       stripchart(bug, vertical=TRUE, axes=FALSE)
			plot( keptCounts ~ factor(sex), ylab = names[index],main = graphMain )
                        points(factor(sex), keptCounts)
#                        plot ( shannon ~ factor(ac) )
#			plot ( keptCounts ~ factor(ac) )
 #                       par(mar=c(1,3,1,2)+0.1)
        plot ( keptCounts ~ factor(ac) )
        points(factor(ac), keptCounts)
#                        par(mar=c(1,3,1,2)+0.1)
        plot( keptCounts ~ factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))))
        points(factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))), keptCounts)
#			plot( keptCounts ~ factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))))

#                        plot( bug ~ factor(experiment) )

#                        plot( bug ~ factor(shannon) )

#                        plot( bug ~ factor(keptCounts) )
#        plot ( keptCounts ~ factor(experiment) )
#        points(factor(experiment), keptCounts)

#                        plot ( bug ~ factor(treatment) )
                                        #                        plot( bug ~ factor(keptCounts), ylab=names[index])
#                        plot( shannon ~ factor(cage), ylab=names[index])
#                         plot( keptCounts ~ factor(cage), ylab=names[index])
#                        par(mar=c(2,3,2,2)+0.1)
        plot( keptCounts ~ factor(cage), ylab=names[index])
        points(factor(cage), keptCounts)
			index=index+1

		}
#	dFrame <- data.frame( names, pValuesSex, pValuesAcuteChronic, pValuesCage, iccCage)#,
	dFrame <- data.frame( names, pValuesSex, pValuesAcuteChronic, pValuesCage, iccCage)#, pValuesExperiment)#, pValuesTreatment)#, pValuesExperiment)
	dFrame <- dFrame [order(pValuesAcuteChronic),]
	dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
	dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
#        dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
#        dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )

	dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
#	dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )
#	dFrame$adjustedShannon <- p.adjust( dFrame$pValuesShannon, method = "BH" )
#        dFrame$adjustedKeptCounts <- p.adjust( dFrame$pValuesKeptCounts, method = "BH" )

#	write.table(dFrame, file=paste("pValuesForTaxa_Sex_cage_acuteChronic_againstShannon", taxa, ".txt",sep=""),
#        print(dim(dFrame))
        write.table(dFrame, file=paste("pValuesForTaxa_keptCounts_Sex_cage_acuteChronic_", taxa, ".txt",sep=""),
#	write.table(dFrame, file=paste("pValuesForTaxa_Sex_cage_acuteChronic_againstkeptCounts", taxa, ".txt",sep=""),
#	write.table(dFrame, file=paste("pValuesForTaxa_Sex_keptCounts_acuteChronic_", taxa, ".txt",sep=""),
        sep="\t",row.names=FALSE)
	dev.off()
}
