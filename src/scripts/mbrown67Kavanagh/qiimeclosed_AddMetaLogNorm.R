## I'm not sure if this is for merged reads or just the forward reads
## Remember that the cutoff is based on the phylum classification...
## Clear and load packages
rm(list=ls())

## Set base directory
baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging"

## Assignment type
dataType <- "closedQIIMER1"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, dataType, "data", sep="/")
processedDir <- paste(baseDir, dataType, "processed", sep="/")
metadataFileName <- "IntAgeMetadata.txt"
baseDataFileName <- "closed_reference_otu_table_L"
## Load metadata
setwd(metadataDir)
metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric")
myMetadata <- read.table(metadataFileName, header=TRUE, sep="\t", colClasses=metaColClasses, skip=0)
## myMetadata <- myMetadata[-1,]
## colnames(myMetadata)[1] <- "MatchFile"
## To deal with conversion of names since they start with numbers
rownames(myMetadata) <- paste0("X", myMetadata[,1])
myMetadata<- myMetadata[,-1]

## Will later be able to go through this and work on the other datasets.
## This code can be reused to get through the qiime taxonomic levels.
taxaLevels <- c(2:7)
for (taxa in taxaLevels){
    ## taxa <- 2
    ## Change over to data directory
    setwd(dataDir)
    inFileName <- paste(baseDataFileName, taxa, ".txt", sep="")
    myT <- read.table(inFileName, sep="\t", header=TRUE, comment.char="@", skip=1)
    numCols <- ncol(myT)
    myColClasses <- c("character", rep("numeric", numCols - 1))
    myT <-read.table(inFileName, header=TRUE, colClasses=myColClasses, sep="\t", comment.char="@", skip=1)

    myTcolnames <- myT[,1]
    myT <- myT[,-1]
    ## pivot
    myT <- as.data.frame(t(myT))
    colnames(myT) <- myTcolnames
    ## colnames(myT) <- as.character(unlist(myT[1,]))
    ## myT <- myT[-1,]
    ## saveRownames <- rownames(myT)
    ## This messes something up...
    ## myT <- as.data.frame(as.matrix(sapply(myT,as.numeric)))
    ## rownames(myT) <- saveRownames

    Sample.Name<-gsub("\\.", "-", unlist(lapply(strsplit(rownames(myT), split="_"),"[[",1)))
    ## Why does the numeric code change?
    Sample.Name<-gsub("X3754", "X3612", Sample.Name)
    rownames(myT) <- Sample.Name

    ## This is where the trouble comes in because they aren't necessary in the same order
    ## counts
    ## myMetadata$counts <- rowSums(myT)
    ## numMetadataCols <- dim(myMetadata)[2]
    nSamples <- dim(myT)[1] ## N total number of samples
    sizeSamples <- rowSums(myT)
    ## number of sequences in a sample
    totalSize <- sum(sizeSamples) ## is total number of counts in the table
    if (taxa == 2){
        sampleDepth <- sizeSamples
    }
    ## Size of each level; it is redundant at the phylum level
    depthAtLevel <- sizeSamples

    myT <- cbind(sampleDepth, depthAtLevel, myT)
    ## combine
    myLogNorm <- myT[,]

    for(i in 1:dim(myT)[1])
    {
        ##Offset for the depth columns
        for(j in 3:dim(myT)[2])
        {
            myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
        }
    }

    myT <- myLogNorm
    ## myT$MatchFile <- as.character(rownames(myT))
    ## myT$MatchFile <- as.factor(myT$MatchFile)
    ## names(myT) <- colnames(myT)
    ## Sample.Name<-gsub("\\.", "-", unlist(lapply(strsplit(rownames(myT), split="_"),"[[",1)))
    ## Why does the numeric code change?
    ## Sample.Name<-gsub("X3754", "X3612", Sample.Name)
    ## Then replace "." with "_"
    ## readDirection <- sapply(strsplit(rownames(myT), split="_"),"[[", 4)
    ## myT<-cbind(myT, Sample.Name)##, readDirection)
    ## myT<-as.data.frame(myT, stringsAsFactors=FALSE)

    merged<-merge(x = myMetadata, y = myT, by = 0, all=TRUE)

    ## Change over to processed data directory
    setwd(processedDir)
    ## output to txt file
    colnames(merged)[1] <- "Sample.ID"
    write.table(merged, file=paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep=""), row.names = FALSE, sep="\t")
}
