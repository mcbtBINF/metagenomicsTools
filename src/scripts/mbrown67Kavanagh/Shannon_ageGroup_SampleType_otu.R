rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
## Below is what I need for the Shannon diversity.
library("vegan")

## for experiment, treatment, and batch, and no confounder

## setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/rdpClassifications")

## taxaLevels <- c( "phylum", "class", "order", "family", "genus")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/otus/")
##taxaLevels <- c( "phylum", "class", "order", "family", "genus")
taxaLevels <- c("otu")
indexS <- 1
ShannonPcage <- list()
ShannonPmodel <- list()
ShannonSummary <- list()

## tissueKept <- c("LI Lumen", "LI Mucosa", "Feces")
## tissueKept <- "LI Lumen"
## tissueKept <- "LI Mucosa"
## tissueKept <- "Feces"
tissueKept <- c("LI Lumen", "LI Mucosa")

for(taxa in taxaLevels )
{

    inFileName <- paste( taxa, "_LogNormwithMetadata_otu.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    numMetadataCols <- 6
    ## Reprocessing for correct interpretation
    ## Somewhat arbitrarily here
    names <- vector()
    pValuesSex<- vector()
    pValuesGroup<- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    iccKeptCounts <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()

    index <- 1
    pdf( paste(taxa, "_testShannonplot_ageGroup_LumenMucosa_otu.pdf", sep=""))
    myT<-myT[myT$Sample.Type %in% tissueKept,]
    animal <- myT$Animal.ID
    sampType <- myT$Sample.Type
    cage <- myT$Pen.Location
    group <- myT$Group
    age <- myT$Age

    myT$Shannon <- apply(myT[,2:(ncol(myT)-numMetadataCols)],1,diversity)
    Shannon <- myT$Shannon
    myFrame <- data.frame( age, group, sampType, cage, animal, Shannon)

    fullModel <- gls( Shannon ~ group + sampType, method="REML",correlation=corCompSymm(form=~1|cage),	data = myFrame )
    reducedModel <- gls( Shannon~ group + sampType, method="REML", data = myFrame )
    fullModelLME <- lme( Shannon~ group + sampType, method="REML", random = ~1|cage, data = myFrame)

    ShannonPcage[indexS] <- anova(fullModel, reducedModel)$"p-value"[2]
    ## Indicates that the cage effect is not significant in explaining the trend in Shannon diveristy
    Shannonlm <- lm(Shannon ~ group + sampType, x=TRUE)
    ShannonPmodel[indexS] <- list(summary(Shannonlm)$coefficients[-1,4])
    plot(group, Shannon, main=taxa, ylab="Shannon Diversity")
    plot(sampType, Shannon, main=taxa, ylab="Shannon Diversity")
    ##list(summary(Shannonlm)$coefficients[,4][-1])
    indexS <- indexS + 1
    ## #        ShannonSummary[indexS]<-summary(Shannonlm)

    ## #        names[indexS] = names(myT)[i]

    ##         indexS <- indexS + 1

    ##         # Building the data.frames to eventually print out the p-values
    ##         DayPatientPV.df <- data.frame(ShannonP)
    ##         DayPatientPV.df <- t(DayPatientPV.df)
    ##         dFrameDayPatient <- DayPatientPV.df


    ## 	graphMain =  paste( names(myT)[i])
    ## 	par(mfrow=c(3,1),
    ##             oma = c(1,1,0,0) + 0.1,
    ##             mar = c(1,4,2.5,0) + 0.1)
    ##         plot(as.numeric(factor(multiWay)), bug, col=as.numeric(batch), pch=16)
    ##         plot(jitter(as.numeric(sl),.3), bug, col=as.numeric(batch), pch=16 + as.numeric(treatment))
    ##         plot(as.numeric(batch), bug, col=as.numeric(sex), pch=16)

    ##         index=index+1
    ## }
    dev.off()
}

dFrame <- data.frame(matrix(unlist(ShannonPmodel), nrow=5, byrow=TRUE))
colnames(dFrame) <- names(ShannonPmodel[[1]])
write.table(dFrame,  file ="pValuesShannon_ageGroup_LumenMucosa_otu.txt", row.names=FALSE, sep="\t")
