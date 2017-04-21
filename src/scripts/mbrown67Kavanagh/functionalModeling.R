## Does the modeling
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
    tissues <- unique(myT$Sample.Type)
    for(tissue in tissues){
        subT <- myT[myT$Sample.Type == tissue,]

        names <- vector()
        pValuesfilePrefix<- vector()
        pValuesfilePrefixWilcox<- vector()
        pValuesCage<- vector()
        iccCage <- vector()
        pValuesSampleType <- vector()
        pValuesAge<- vector()
        pValuesAnimal<- vector()
        pValuesInteraction <- vector()
        allpvals <- list()

        index <- 1

        ## Over the non-metadata columns
        ## Will replace with a function call without the loop
        ## Add in some timing information.
        for ( j in c((endMetadataIndex + 1):ncol(subT))) {
            bug <- subT[,j]
            ## Removing rare taxa
            ## This should be a bit more flexible
            if( sum(bug != 0 ) > nrow(subT) / 4 ) {
                names[index] <- names(subT)[j]

                varCare <- factor(unlist(subT[filePrefix]))
                sampleType <- factor(subT$Sample.Type)

                cage <-  factor(subT$Pen.Location)

                myFrame <- data.frame(bug, varCare, sampleType, cage)

                fullModel <- lm(bug~varCare, data= myFrame)
                mixedAnova <- anova(fullModel)
                ## This method of getting p-values only works for the univariate case
                pValuesfilePrefix[index] <- mixedAnova$"Pr(>F)"[1]
                pValuesfilePrefixWilcox[index] <- wilcox.test( bug[subT[filePrefix] == "Young"],    bug[subT[filePrefix] == "Old"])$p.value

                index <- index + 1
            }
        }

        dFrame <- data.frame(names, pValuesfilePrefix, pValuesfilePrefixWilcox)
        rownames(dFrame) <- dFrame$names
        dFrame <- dFrame [order(dFrame$pValuesfilePrefix),]
        dFrame$adjustedpValuesfilePrefix <- p.adjust( dFrame$pValuesfilePrefix, method = "BH" )
        dFrame$adjustedpValuesfilePrefixWilcox <- p.adjust( dFrame$pValuesfilePrefixWilcox, method = "BH" )
        oldColnames <- colnames(dFrame)
        ## This can probably be better written
        colnames(dFrame) <- c("names", paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"))
        write.table(dFrame, file=paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, ".txt",sep=""), sep="\t", row.names = FALSE)
        colnames(dFrame) <- oldColnames
        pdf( paste(dataType, "_L_", taxa, "_", tissue, "_", filePrefix, "_boxplots.pdf", sep=""))
        index <- 1
        for (nameIter in rownames(dFrame)) {
            ## This should be stored once near the top and altered there
            par(mfrow=c(1,1),
                oma = c(0, 3, 2, 0) + 0.1,
                mar = c(3, 2, 3, 1) + 0.1,
                mgp = c(3, 1.25, 0))

            bug <- subT[,nameIter]
            varCare <- factor(unlist(subT[filePrefix]))
            sampleType <- factor(subT$Sample.Type)
            cage <-  factor(subT$Pen.Location)
            myFrame <- data.frame(bug, varCare, sampleType, cage)

            boxplot( bug ~ varCare , ylab="Log normalized abundance", pch=16,
                    main = paste(nameIter, filePrefix,"\n p-value linear model", format(dFrame$adjustedpValuesfilePrefix[index],digits=3), "\n p-value Wilcox test", format(dFrame$adjustedpValuesfilePrefixWilcox[index],digits=3) ) )

            stripchart(bug ~ varCare ,
                       data = myFrame, vertical = TRUE, pch = 16, add=TRUE, ylab = "Log normalized abundance", method="jitter")

            index = index + 1
        }

        hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for\n", filePrefix, " at Taxonomic Level:", taxa),
             ylab = "Count", xlab="Unadjusted p-values")
        hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for\n", filePrefix, " at Taxonomic Level:", taxa),
             ylab = "Count", xlab="Unadjusted p-values")

        dev.off()
    }
}

##for (taxa in taxaLevels){
for(tissue in tissues){
    pdf( paste(dataType, "_", tissue, "_", filePrefix, "_ShannonDiversity.pdf", sep=""))
    par(mfrow=c(1,1),
        oma = c(0, 3, 2, 0) + 0.1,
        mar = c(3, 2, 3, 1) + 0.1,
        mgp = c(3, 1.25, 0))


    names <- vector()
    pValuesfilePrefix<- vector()
    pValuesfilePrefixWilcox <- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    pValuesSampleType <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()
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

        subT <- myT[myT$Sample.Type == tissue,]

        subT$ShannonDiversity <- apply(subT[,(endMetadataIndex+1):ncol(subT)],1,diversity)
        bug <- subT$ShannonDiversity
        ShannonDiversityList[[index]] <- bug

        varCare <- factor(unlist(subT[filePrefix]))
        sampleType <- factor(subT$Sample.Type)
        cage <-  factor(subT$Pen.Location)
        myFrame <- data.frame(bug, varCare, sampleType, cage)

        fullModel <- lm(bug~varCare, data= myFrame)
        mixedAnova <- anova(fullModel)
        pValuesfilePrefix[index] <- mixedAnova$"Pr(>F)"[1]
        pValuesfilePrefixWilcox[index] <- wilcox.test( bug[subT[filePrefix] == "Young"],    bug[subT[filePrefix] == "Old"])$p.value

        boxplot( bug ~ varCare , ylab="Shannon Diversity", pch=16,
                main = paste(filePrefix, " linear model p-value", format(pValuesfilePrefix[index],digits=3),"\n Wilcox test p-value", format(pValuesfilePrefixWilcox[index],digits=3), "\n at taxonomic level ", taxa ) )

        stripchart(bug ~ varCare ,
                   data = myFrame,vertical = TRUE, pch = 16, add=TRUE, ylab = "Shannon Diversity", method="jitter")

        index <- index + 1
    }

    ## Will have to do something about labeling the taxonomic levels
    ## hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for ", filePrefix, "\nat Taxonomic Level:", taxa))
    ## hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for ", filePrefix, "\nat Taxonomic Level:", taxa))

    dev.off()

    dSFrame <- data.frame(pValuesfilePrefix, pValuesfilePrefixWilcox)
    dSFrame <- dSFrame[order(dSFrame$pValuesfilePrefix),]
    dSFrame$adjustedpValuesfilePrefix <- p.adjust( dSFrame$pValuesfilePrefix, method = "BH" )
    dSFrame$adjustedpValuesfilePrefixWilcox <- p.adjust(dSFrame$pValuesfilePrefixWilcox, method = "BH")
    rownames(dSFrame) <- dSFrame$names
    colnames(dFrame) <- c(paste0("p-values",filePrefix), paste0("p-values",filePrefix,"Wilcox"), paste0("adjusted_p-values",filePrefix), paste0("adjustedp-values",filePrefix,"Wilcox"))

    write.table(dSFrame, file=paste(dataType, tissue, "_", filePrefix, "_ShannonDiversity.txt", sep=""), sep="\t", row.names=FALSE)
}
## There's really no need to go through this trouble as Shannon diversity is not usually MHC'd
## pdf( paste(dataType, "_", tissue, "_", filePrefix, "_ShannonDiversity.pdf", sep=""))
## index <- 1
## for (nameIter in rownames(dSFrame)) {
##     ## Handle this at the top
##     par(mfrow=c(1,1),
##         oma = c(0, 3, 2, 0) + 0.1,
##         mar = c(3, 2, 2, 1) + 0.1,
##         mgp = c(3, 1.25, 0))

##     bug <- ShannonDiversityList[[index]]

##     varCare <- factor(unlist(subT[filePrefix]))
##     sampleType <- factor(subT$Sample.Type)

##     cage <-  factor(subT$Pen.Location)

##     myFrame <- data.frame(bug, varCare, sampleType, cage)

##     boxplot( bug ~ varCare , ylab="Shannon Diversity",
##             main = paste("Age filePrefix p-value", format(dSFrame$adjustedpValuesfilePrefix[index],digits=3) ) )

##     stripchart(bug ~ varCare ,
##                data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = "Shannon Diversity")

##     index = index + 1
## }

## hist(pValuesfilePrefix, breaks=20, main=paste0("Unadjusted p-values (parametric linear model) for ", filePrefix, " at Taxonomic Level:", taxa))
## hist(pValuesfilePrefixWilcox, breaks=20, main=paste0("Unadjusted p-values (non-parametric Wilcox test) for ", filePrefix, " at Taxonomic Level:", taxa))

## dev.off()
## }
