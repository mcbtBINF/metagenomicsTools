## Does the modeling on abundances, Shannon diversity and MDS axes
## parametric and nonparametric
## And the Shannon diversity
## And modeling on MDS axes
## Will this be easy to modify with depth?
## Rare taxa filter
## Sequencing Depth as a confounder should be checked out automatically too.

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")
baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging"
setwd(baseDir)

dataType <- "closedQIIMER1"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, dataType, sep="/")
processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(baseDir, dataType, "analysis", sep="/")
baseDataFileName <- "closed_reference_otu_table_L"
## It would be great to make the name change the variables in place
## This should replace the each modeling code
## This also include SAMPLETYPE
filePrefix <- "Group"

taxaLevels <- c(2:7)
for (taxa in taxaLevels){
    setwd(processedDir)
    inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
    myT <- read.table(inFileName, header=TRUE, sep="\t")
    numCols <- ncol(myT)

    ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    ## The 9 could be done better...
    myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
    myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)

    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    setwd(analysisDir)
    ## Need to find a way to do all tissues
    ## No, that is another scripts job as it involves a second variable.

    names <- vector()
    pValuesfilePrefix <- vector()
    pValuesfilePrefixWilcox <- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    iccAnimal <- vector()
    pValuesSampleType <- vector()
    pValuesfilePrefixMixed <- vector()
    pValuesSampleTypeMixed <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()
    ## Additional checks
    pValuesDepth <- vector()
    ## pValuesDepthWilcox <- vector()
    pValuesDepthwithfilePrefix <- vector()
    pValuesfilePrefixwithDepth <- vector()
    pValuesfilePrefixwithSampleType <- vector()

    index <- 1

    ## Over the non-metadata columns
    ## Will replace with a function call without the loop
    ## Add in some timing information.
    for ( j in c((endMetadataIndex + 1):ncol(myT))) {
        bug <- myT[,j]
        ## Removing rare taxa
        ## This should be a bit more flexible
        if( sum(bug != 0 ) > nrow(myT) / 4 ) {
            names[index] <- names(myT)[j]

            varCare <- factor(unlist(myT[filePrefix]))
            sampleType <- factor(myT$Sample.Type)

            cage <-  factor(myT$Pen.Location)
            depth <- myT$depthAtLevel
            animal <- myT$Animal.ID

            myFrame <- data.frame(bug, varCare, sampleType, depth, animal)

            varCareModel <- lm(bug~varCare, data= myFrame)
            depthModel <- lm(bug~depth, data = myFrame)
            depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
            myAnova <- anova(varCareModel)

            varCareModelwSampleType <- lm(bug~varCare + sampleType, data= myFrame)
            sampleTypeModel <- lm(bug~sampleType, data=myFrame)
            fullModel <-
                gls( bug~  varCare + sampleType,
                    method="REML", correlation=corCompSymm(form=~1|factor(animal)),data = myFrame )

            mixedAnova <- anova(fullModel)
            pValuesfilePrefixMixed[index] <- mixedAnova$"p-value"[2]
            pValuesSampleTypeMixed[index] <- mixedAnova$"p-value"[3]

            pValuesAnimal[index] <- anova(lm( bug ~ animal ))$"Pr(>F)"[1]
            iccAnimal[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


            ## This method of getting p-values only works for the univariate case
            pValuesfilePrefix[index] <- myAnova$"Pr(>F)"[1]
            pValuesfilePrefixWilcox[index] <- wilcox.test( bug[myT[filePrefix] == "Young"],    bug[myT[filePrefix] == "Old"])$p.value
            pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
            ## pValuesDepthwithfilePrefix[index] <-
            pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
            pValuesSampleType[index] <- anova(sampleTypeModel)$"Pr(>F)"[1]
            pValuesfilePrefixwithSampleType[index] <- anova(varCareModelwSampleType, sampleTypeModel)$"Pr(>F)"[2]

            ## Model for depth
            ## Model for varCare and depth

            index <- index + 1
        }
    }

    dFrame <- data.frame(names, pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth, pValuesfilePrefixMixed, pValuesSampleTypeMixed, pValuesAnimal, iccAnimal)
    dFrame <- dFrame [order(dFrame$pValuesfilePrefix),]
    rownames(dFrame) <- dFrame$names
    dFrame$adjustedpValuesfilePrefix <- p.adjust( dFrame$pValuesfilePrefix, method = "BH" )
    dFrame$adjustedpValuesfilePrefixWilcox <- p.adjust( dFrame$pValuesfilePrefixWilcox, method = "BH" )
    dFrame$adjustedpValuesDepth <- p.adjust(dFrame$pValuesDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dFrame$pValuesfilePrefixwithDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixMixed <- p.adjust(dFrame$pValuesfilePrefixMixed, method="BH")
    dFrame$adjustedpValuesSampleTypeMixed <- p.adjust(dFrame$pValuesSampleTypeMixed, method="BH")
    dFrame$adjustedpValuesAnimal <- p.adjust(dFrame$pValuesAnimal, method="BH")
    oldColnames <- colnames(dFrame)
    ## This can probably be better written
    colnames(dFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("p-valuesfromMixed", filePrefix), paste0("p-valuesSampleTypefromMixed"), paste0("p-valuesAnimal"), "iccAnimal", paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix), paste0("adjusted_p-valuesfromMixed",filePrefix), paste0("adjusted_p-valuesSampleTypefromMixed"),paste0("adjusted_p-valuesAnimal"))
    write.table(dFrame, file=paste(dataType, "_L_", taxa, "_", filePrefix, "_SampleType_AnimalMixed.txt",sep=""), sep="\t", row.names = FALSE)
    colnames(dFrame) <- oldColnames
    pdf( paste(dataType, "_L_", taxa, "_", filePrefix, "_SampleType_AnimalMixedboxplots.pdf", sep=""))

    index <- 1
    for (nameIter in rownames(dFrame)) {
        ## This should be stored once near the top and altered there
        par(mfrow=c(3,1),
            oma = c(0, 3, 2, 0) + 0.1,
            mar = c(3, 2, 3, 1) + 0.1,
            mgp = c(3, 1.25, 0))

        bug <- myT[,nameIter]
        varCare <- factor(unlist(myT[filePrefix]))
        depth <- myT$depthAtLevel
        sampleType <- factor(myT$Sample.Type)
        cage <-  factor(myT$Pen.Location)
        animal <- factor(myT$Animal.ID)
        myFrame <- data.frame(bug, varCare, sampleType, animal, depth)

        boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=16,
                main = paste(filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

        stripchart(bug ~ varCare ,
                   data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

        plot( bug ~ depth , ylab="Log normalized abundance", pch=16,
             main = paste("Depth p-value linear model", format(dFrame$adjustedpValuesDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

        stripchart(bug ~ depth ,
                   data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

        plot( bug ~ sampleType , ylab="Log normalized abundance", pch=16,
             main = paste(filePrefix, "Sample Type p-value linear model", format(dFrame$adjustedpValuesSampleType[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

        stripchart(bug ~ sampleType ,
                   data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

        ## boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=16,
        ##         main = paste(filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

        ## stripchart(bug ~ varCare ,
        ##            data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")


        boxplot( bug ~ animal, ylab="Log normalized abundance",
                main=paste("Animal p-value", format(dFrame$adjustedpValuesAnimal[index],digits=3)))

        stripchart(bug ~ animal,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])


        mtext(nameIter, outer=TRUE, cex = 0.7)
        mtext("Log normalized abundance", outer=TRUE, side=2, line=1)

        index = index + 1
    }

    hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixwithDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " after depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesSampleType, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "SampleType", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesAnimal, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Animal", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")


    dev.off()
}

## This is for the MDS analysis
for (taxa in taxaLevels){
    setwd(processedDir)
    inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
    myT <- read.table(inFileName, header=TRUE, sep="\t")
    numCols <- ncol(myT)

    ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    ## The 9 could be done better...
    myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
    myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)

    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    setwd(analysisDir)
    ## Need to find a way to do all tissues
    ## No, that is another scripts job as it involves a second variable.

    ## This will need to be replaced
    inFileName <- paste(baseDataFileName, taxa, "_pcoa.txt",sep="")

    myPCoA <- read.table(inFileName, sep="\t", header=TRUE)

    names <- vector()
    pValuesfilePrefix <- vector()
    pValuesfilePrefixWilcox <- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    iccAnimal <- vector()
    pValuesSampleType <- vector()
    pValuesfilePrefixMixed <- vector()
    pValuesSampleTypeMixed <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()
    ## Additional checks
    pValuesDepth <- vector()
    ## pValuesDepthWilcox <- vector()
    pValuesDepthwithfilePrefix <- vector()
    pValuesfilePrefixwithDepth <- vector()
    pValuesfilePrefixwithSampleType <- vector()



    index <- 1

    ## Over the non-metadata columns
    ## Will replace with a function call without the loop
    ## Add in some timing information.
    for ( j in 1:10) {
        bug <- myPCoA[,j]
        ## Removing rare taxa
        ## This should be a bit more flexible
        names[index] <- names(myPCoA)[j]

        varCare <- factor(unlist(myT[filePrefix]))
        sampleType <- factor(myT$Sample.Type)

        cage <-  factor(myT$Pen.Location)
        depth <- myT$depthAtLevel
        animal <- myT$Animal.ID

        myFrame <- data.frame(bug, varCare, sampleType, depth, animal)

        varCareModel <- lm(bug~varCare, data= myFrame)
        depthModel <- lm(bug~depth, data = myFrame)
        depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
        myAnova <- anova(varCareModel)

        varCareModelwSampleType <- lm(bug~varCare + sampleType, data= myFrame)
        sampleTypeModel <- lm(bug~sampleType, data=myFrame)
        fullModel <-
            gls( bug~  varCare + sampleType,
                method="REML", correlation=corCompSymm(form=~1|factor(animal)),data = myFrame )

        mixedAnova <- anova(fullModel)
        pValuesfilePrefixMixed[index] <- mixedAnova$"p-value"[2]
        pValuesSampleTypeMixed[index] <- mixedAnova$"p-value"[3]

        pValuesAnimal[index] <- anova(lm( bug ~ animal ))$"Pr(>F)"[1]
        iccAnimal[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


        ## This method of getting p-values only works for the univariate case
        pValuesfilePrefix[index] <- myAnova$"Pr(>F)"[1]
        pValuesfilePrefixWilcox[index] <- wilcox.test( bug[myT[filePrefix] == "Young"],    bug[myT[filePrefix] == "Old"])$p.value
        pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
        ## pValuesDepthwithfilePrefix[index] <-
        pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
        pValuesSampleType[index] <- anova(sampleTypeModel)$"Pr(>F)"[1]
        pValuesfilePrefixwithSampleType[index] <- anova(varCareModelwSampleType, sampleTypeModel)$"Pr(>F)"[2]

        index <- index + 1
    }

    dFrame <- data.frame(names, pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth, pValuesfilePrefixMixed, pValuesSampleTypeMixed, pValuesAnimal, iccAnimal)
    dFrame <- dFrame [order(dFrame$pValuesfilePrefix),]
    rownames(dFrame) <- dFrame$names
    dFrame$adjustedpValuesfilePrefix <- p.adjust( dFrame$pValuesfilePrefix, method = "BH" )
    dFrame$adjustedpValuesfilePrefixWilcox <- p.adjust( dFrame$pValuesfilePrefixWilcox, method = "BH" )
    dFrame$adjustedpValuesDepth <- p.adjust(dFrame$pValuesDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dFrame$pValuesfilePrefixwithDepth, method="BH")
    dFrame$adjustedpValuesfilePrefixMixed <- p.adjust(dFrame$pValuesfilePrefixMixed, method="BH")
    dFrame$adjustedpValuesSampleTypeMixed <- p.adjust(dFrame$pValuesSampleTypeMixed, method="BH")
    dFrame$adjustedpValuesAnimal <- p.adjust(dFrame$pValuesAnimal, method="BH")
    oldColnames <- colnames(dFrame)
    ## This can probably be better written
    colnames(dFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("p-valuesfromMixed", filePrefix), paste0("p-valuesSampleTypefromMixed"), paste0("p-valuesAnimal"), "iccAnimal", paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix), paste0("adjusted_p-valuesfromMixed",filePrefix), paste0("adjusted_p-valuesSampleTypefromMixed"),paste0("adjusted_p-valuesAnimal"))

    write.table(dFrame, file=paste(dataType, "_L_", taxa, "_", filePrefix, "_SampleType_Animal_MixedMDS.txt",sep=""), sep="\t", row.names = FALSE)
    colnames(dFrame) <- oldColnames
    pdf( paste(dataType, "_L_", taxa, "_", filePrefix, "_SampleType_Animal_MixedMDSboxplots.pdf", sep=""))
    index <- 1
    for (j in 1:10) {
        ## This should be stored once near the top and altered there
        par(mfrow=c(3,1),
            oma = c(0, 3, 2, 0) + 0.1,
            mar = c(3, 2, 3, 1) + 0.1,
            mgp = c(3, 1.25, 0))

        bug <- myPCoA[,j]
        varCare <- factor(unlist(myT[filePrefix]))
        depth <- myT$depthAtLevel
        sampleType <- factor(myT$Sample.Type)
        cage <-  factor(myT$Pen.Location)
        myFrame <- data.frame(bug, varCare, sampleType, cage, depth)

        boxplot( bug ~ varCare , ylab="MDS Axis", pch=16,
                main = paste(filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

        stripchart(bug ~ varCare ,
                   data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "MDS Axis", method="jitter")

        plot( bug ~ depth , ylab="MDS Axis", pch=16,
             main = paste("Depth p-value linear model", format(dFrame$adjustedpValuesDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

        plot( bug ~ varCare + depth , ylab="MDS Axis", pch=16,
             main = paste(filePrefix, "and Depth p-value linear model", format(dFrame$adjustedpValuesfilePrefixwithDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

        plot( bug ~ sampleType , ylab="Log normalized abundance", pch=16,
             main = paste(filePrefix, "Sample Type p-value linear model", format(dFrame$adjustedpValuesSampleType[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

        stripchart(bug ~ sampleType ,
                   data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

        ## boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=16,
        ##         main = paste(filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

        ## stripchart(bug ~ varCare ,
        ##            data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")


        boxplot( bug ~ animal, ylab="Log normalized abundance",
                main=paste("Animal p-value", format(dFrame$adjustedpValuesAnimal[index],digits=3)))

        stripchart(bug ~ animal,
                   data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])


        ## mtext(nameIter, outer=TRUE, cex = 0.7)
        mtext(paste0("MDS Axis ", j), outer=TRUE, side=2, line=1)

        ## stripchart(bug ~ varCare + depth ,
        ##           data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "MDS Axis", method="jitter")
        ## mtext(nameIter, outer=TRUE, cex = 0.7)
        ## mtext(paste0("MDS Axis ", j), outer=TRUE, side=2, line=1)

        index = index + 1
    }

    hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for\n", filePrefix, " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesfilePrefixwithDepth, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " after depth", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesSampleType, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "SampleType", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")
    hist(pValuesAnimal, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", "Animal", " at Taxonomic Level:", taxa),
         ylab = "Count", xlab="Unadjusted p-values")

    dev.off()
}

##for (taxa in taxaLevels){
pdf( paste(dataType, "_", filePrefix, "_SampleType_Animal_MixedShannonDiversity.pdf", sep=""))
par(mfrow=c(3,1),
    oma = c(0, 3, 2, 0) + 0.1,
    mar = c(3, 2, 3, 1) + 0.1,
    mgp = c(3, 1.25, 0))

names <- vector()
pValuesfilePrefix <- vector()
pValuesfilePrefixWilcox <- vector()
pValuesCage<- vector()
iccCage <- vector()
iccAnimal <- vector()
pValuesSampleType <- vector()
pValuesfilePrefixMixed <- vector()
pValuesSampleTypeMixed <- vector()
pValuesAge<- vector()
pValuesAnimal<- vector()
pValuesInteraction <- vector()
allpvals <- list()
## Additional checks
pValuesDepth <- vector()
## pValuesDepthWilcox <- vector()
pValuesDepthwithfilePrefix <- vector()
pValuesfilePrefixwithDepth <- vector()
pValuesfilePrefixwithSampleType <- vector()

ShannonDiversityList <- list()

index <- 1
for(taxa in taxaLevels){
    setwd(processedDir)
    inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
    myT <- read.table(inFileName, header=TRUE, sep="\t")
    numCols <- ncol(myT)

    ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    ## The 9 could be done better...
    myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
    myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)

    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    setwd(analysisDir)

    myT$ShannonDiversity <- apply(myT[,(endMetadataIndex+1):ncol(myT)],1,diversity)
    bug <- myT$ShannonDiversity
    ShannonDiversityList[[index]] <- bug

    varCare <- factor(unlist(myT[filePrefix]))
    sampleType <- factor(myT$Sample.Type)
    cage <-  factor(myT$Pen.Location)
    animal <- factor(myT$Animal.ID)
    depth <- myT$depthAtLevel
    myFrame <- data.frame(bug, varCare, sampleType, cage, depth, animal)
    ## Insert p-value assignments here
    varCareModel <- lm(bug~varCare, data= myFrame)
    depthModel <- lm(bug~depth, data = myFrame)
    depthWithCareModel <- lm(bug~varCare + depth, data=myFrame)
    myAnova <- anova(varCareModel)

    varCareModelwSampleType <- lm(bug~varCare + sampleType, data= myFrame)
    sampleTypeModel <- lm(bug~sampleType, data=myFrame)
    fullModel <-
        gls( bug~  varCare + sampleType,
            method="REML", correlation=corCompSymm(form=~1|factor(animal)),data = myFrame )

    mixedAnova <- anova(fullModel)
    pValuesfilePrefixMixed[index] <- mixedAnova$"p-value"[2]
    pValuesSampleTypeMixed[index] <- mixedAnova$"p-value"[3]

    pValuesAnimal[index] <- anova(lm( bug ~ animal ))$"Pr(>F)"[1]
    iccAnimal[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


    ## This method of getting p-values only works for the univariate case
    pValuesfilePrefix[index] <- myAnova$"Pr(>F)"[1]
    pValuesfilePrefixWilcox[index] <- wilcox.test( bug[myT[filePrefix] == "Young"],    bug[myT[filePrefix] == "Old"])$p.value
    pValuesDepth[index] <- anova(depthModel)$"Pr(>F)"[1]
    ## pValuesDepthwithfilePrefix[index] <-
    pValuesfilePrefixwithDepth[index] <- anova(depthWithCareModel, depthModel)$"Pr(>F)"[2]
    pValuesSampleType[index] <- anova(sampleTypeModel)$"Pr(>F)"[1]
    pValuesfilePrefixwithSampleType[index] <- anova(varCareModelwSampleType, sampleTypeModel)$"Pr(>F)"[2]

    boxplot( bug ~ varCare , ylab="Shannon Diversity", pch=16,
            main = paste(filePrefix, " linear model p-value", format(pValuesfilePrefix[index],digits=3),"\n Wilcox test p-value", format(pValuesfilePrefixWilcox[index],digits=3), "\n at taxonomic level ", taxa ) )

    stripchart(bug ~ varCare ,
               data = myFrame,vertical = TRUE, pch = 16, add=TRUE, ylab = "Shannon Diversity", method="jitter")
    plot( bug ~ depth , ylab="Shannon Diversity", pch=16,
         main = paste("Depth p-value linear model", format(dSFrame$adjustedpValuesDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

    stripchart(bug ~ depth ,
               data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

    plot( bug ~ varCare + depth , ylab="Shannon diversity", pch=16,
         main = paste(filePrefix, "and Depth p-value linear model", format(dSFrame$adjustedpValuesfilePrefixwithDepth[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

    stripchart(bug ~ varCare + depth ,
               data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Shannon Diversity", method="jitter")

    plot( bug ~ sampleType , ylab="Log normalized abundance", pch=16,
         main = paste(filePrefix, "Sample Type p-value linear model", format(dSFrame$adjustedpValuesSampleType[index],digits=3)), col=ifelse(varCare == "Young", "blue", "red"))

    stripchart(bug ~ sampleType ,
               data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

    ## boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=16,
    ##         main = paste(filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

    ## stripchart(bug ~ varCare ,
    ##            data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")


    boxplot( bug ~ animal, ylab="Log normalized abundance",
            main=paste("Animal p-value", format(dSFrame$adjustedpValuesAnimal[index],digits=3)))

    stripchart(bug ~ animal,
               data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

    mtext(paste0("Shannon diversity at level ",taxa), outer=TRUE, cex = 0.7)
    mtext("Shannon Diversity", outer=TRUE, side=2, line=1)


    index <- index + 1
}

## Will have to do something about labeling the taxonomic levels
## hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for ", filePrefix, "\nat Taxonomic Level:", taxa))
## hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for ", filePrefix, "\nat Taxonomic Level:", taxa))

dev.off()

dSFrame <- data.frame(taxaLevels, pValuesfilePrefix, pValuesfilePrefixWilcox, pValuesDepth, pValuesfilePrefixwithDepth, pValuesfilePrefixMixed, pValuesSampleTypeMixed, pValuesAnimal, iccAnimal)

dSFrame <- dSFrame[order(dSFrame$pValuesfilePrefix),]
dSFrame$adjustedpValuesfilePrefix <- p.adjust( dSFrame$pValuesfilePrefix, method = "BH" )
dSFrame$adjustedpValuesfilePrefixWilcox <- p.adjust(dSFrame$pValuesfilePrefixWilcox, method = "BH")
rownames(dSFrame) <- dSFrame$names
dSFrame$adjustedpValuesDepth <- p.adjust(dSFrame$pValuesDepth, method="BH")
dSFrame$adjustedpValuesfilePrefixwithDepth <- p.adjust(dSFrame$pValuesfilePrefixwithDepth, method="BH")
dSFrame$adjustedpValuesfilePrefixMixed <- p.adjust(dSFrame$pValuesfilePrefixMixed, method="BH")
dSFrame$adjustedpValuesSampleTypeMixed <- p.adjust(dSFrame$pValuesSampleTypeMixed, method="BH")
dSFrame$adjustedpValuesAnimal <- p.adjust(dSFrame$pValuesAnimal, method="BH")

oldColnames <- colnames(dSFrame)
## This can probably be better written

colnames(dSFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("p-valuesDepth"), paste0("p-valuesDepthwith", filePrefix), paste0("p-valuesfromMixed", filePrefix), paste0("p-valuesSampleTypefromMixed"), paste0("p-valuesAnimal"), "iccAnimal", paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"), paste0("adjusted_p-valuesDepth"), paste0("adjusted_p-valuesDepthwith", filePrefix), paste0("adjusted_p-valuesfromMixed",filePrefix), paste0("adjusted_p-valuesSampleTypefromMixed"),paste0("adjusted_p-valuesAnimal"))

write.table(dSFrame, file=paste(dataType, "_", filePrefix, "_SampleType_Animal_Mixed_ShannonDiversity.txt", sep=""), sep="\t", row.names=FALSE)

