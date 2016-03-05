rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")

# for experiment, treatment, and batch, and no confounder

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/FodorProcessedForwardReadsNoCutoff/noLowTaxas")
#setwd("C:\\MarkRatDataAug2015")

taxaLevels <- c( "p", "c", "o", "f", "g", "otu");
# All of this data should be present in any analysis, so there is no reason to make it optional or subject to a switch/case statement.
raw.Freads<-read.table("../../qiime18_Freads_cr.txt", header=TRUE, comment.char="@", sep="\t")
colcounts<-colSums(raw.Freads[,2:70])
#boxplot(colcounts)
#points(rnorm(69,mean=1,sd=0.05), colcounts)
# Drop low reads (less than 1000) and negative controls and water and technical replicates
excludedCols<-c("X.OTUID", "X45", "X50", "X59tr", "X60tr", "X80", "X81", "X82", "X83", "Water", "N1", "N2", "N3", "P2", "taxonomy",
# The next lines get rid of batch03
#                "X70", "X71", "X72", "X73", "X74", "X75", "X76", "X77", "X78", "X79", "X80", "X81", "X82", "X83", "X84", "X85", "X86", "X87", "X88", "X89", "X90","X91",
                )
# Can drop AM to compare CF to AF
# "X30", "X31",  "X32",  "X33",  "X34",  "X35",  "X36",  "X37",  "X38",  "X39",  "X40",  "X41",  "X42", "X43", "X44")
# Can drop CF to compare AM to AF
# "X46",  "X47",  "X48",  "X49",  "X50",  "X51",  "X52",  "X52",  "X53",  "X54",  "X55",  "X56",  "X57",  "X58",  "X59",  "X60",  "X61")

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
        pcoaFile <- paste("pcoa_",taxa,".txt",sep="")
        myMDS <-read.table(pcoaFile,header=TRUE,sep="\t")

	# squash the negative controls
	myT <- myT[ which( myT$sex != "Non"), ]
        myT$shannon<-diversity(inforShannon[,8:(numColsS - 1)], index="shannon", MARGIN=1)
#        myT$keptCounts<-keptCounts

	# our initial model not worrying about confounders except cage
	names <- vector()
	pValuesSex<- vector()
	pValuesAcuteChronic<- vector()
	pValuesCage<- vector()
        iccCage <- vector()
        iccKeptCounts <- vector()
        pValuesBatch<- vector()
        pValuesTreatment<- vector()
        pValuesShannon<- vector()
        pValuesKeptCounts<- vector()
	# WARNING:  EXPERIMENT IS CONFOUNDED WITH AC  + SEX - INTERPRET WITH CAUTION!!!!!
	pValuesExperiment <- vector()
	index <- 1
	pdf( paste(taxa, "boxplots.pdf", sep=""))

	for( i in 1:20) #9:ncol(myT))
            #This is the poor OTU/resolution limiting step
                                        #            for()
            #Still does this based on the right count of bugs
            if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
                    {
                                        #                        for( j in 1:4){
                        #reassigns bug to MDS value
                            bug <- myMDS[,i]#myT[,i]
                            ac <- myT$acuteOrChronic
                            sex <- myT$sex
                            cage <- myT$cage
                            experiment <- myT$expriment
                            batch <- myT$batch
                            treatment <- myT$treatment
                            shannon <- myT$shannon
                            keptCounts <- myT$keptCounts

                            myFrame <- data.frame(bug, ac, sex, cage, experiment, batch, treatment)


                            fullModel <- gls( bug~  ac + sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),				data = myFrame )
                            reducedModel <- gls( bug~  ac + sex, method="REML",	data = myFrame )
                            fullModelLME <- lme(bug~  ac + sex, method="REML", random = ~1|factor(cage), data = myFrame)

# It seems like reducing the anova calls could speed things up a bit
                            pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[2]
                            pValuesSex[index] <- anova(fullModelLME)$"p-value"[3]
                            pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[5]
                            pValuesBatch[index] <- anova(fullModelLME)$"p-value"[5]
                            pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[5]
                            pValuesShannon[index] <- anova(fullModelLME)$"p-value"[5]
                                        #Why do you mix model functions here?
                            pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
                            iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

                            names[index] = names(myMDS)[i]

                            graphMain =  paste( names(myMDS)[i], "\n",
                                " pSex=", format( pValuesSex[index], digits=3),
                                " pAcuteChronic= ", format( pValuesAcuteChronic[index],digits=3),
#                            " pBatch= ", format(pValuesBatch[index], digits=3), "\n",
#                            " pTreatment= ", format(pValuesTreatment[index], digits=3),
#                            " pExperiment= ", format(pValuesExperiment[index], digits=3),
                            #" pShannon= ", format(pValuesShannon[index], digits=3),
                            #" pKeptCounts= ", format(pValuesKeptCounts[index], digits=3),
                                " pCage= " , format( pValuesCage[index], digits=3),
                                " icc= " , format( iccCage[index], digits=3 ), sep="")
                            par(mfrow=c(4,1),
                                oma = c(1,1,0,0) + 0.1,
                                mar = c(1,4,2.5,0) + 0.1)

                            plot( bug ~ factor(sex), ylab = names[index],main = graphMain )
                            points(factor(sex), bug)

                            plot ( bug ~ factor(ac) )
                            points(factor(ac), bug)

                            plot( bug ~ factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))))
                            points(factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))), bug)

#                        plot ( bug ~ factor(experiment) )
#                        points(factor(experiment), bug)

                            plot( bug ~ factor(cage), ylab=names[index])
                            points(factor(cage), bug)
                            index=index+1

		}

	dFrame <- data.frame( names, pValuesSex, pValuesAcuteChronic, pValuesCage, iccCage)#, pValuesExperiment)#, pValuesBatch)#, pValuesExperiment)
	dFrame <- dFrame [order(pValuesAcuteChronic),]
	dFrame$adjustedAcuteChronic <- p.adjust( dFrame$pValuesAcuteChronic, method = "BH" )
	dFrame$adjustedSex<- p.adjust( dFrame$pValuesSex, method = "BH" )
#        dFrame$adjustedBatch<- p.adjust( dFrame$pValuesBatch, method = "BH" )
#        dFrame$adjustedTreatment<- p.adjust( dFrame$pValuesTreatment, method = "BH" )

	dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
#	dFrame$adjustedExperiment <- p.adjust( dFrame$pValuesExperiment, method = "BH" )

        write.table(dFrame, file=paste("pValuesForTaxa_MDS1to20_Sex_cage_acuteChronic_", taxa, ".txt",sep=""),
        sep="\t",row.names=FALSE)
	dev.off()
}
